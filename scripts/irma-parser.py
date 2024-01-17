"""
IRMA Helpers - parse output directory
"""
import json
import subprocess
from pathlib import Path
from collections import defaultdict

import typer
import toyplot
import pandas as pd
import numpy as np
from toyplot import svg, pdf, png, html


app = typer.Typer(
    help="irmakit - helper functions for parsing irma output",
    add_completion=False,
    pretty_exceptions_short=True,
    no_args_is_help=True,
    rich_markup_mode="rich",
)


@app.command(no_args_is_help=True)
def parse(
    irmadir: Path = typer.Option(
        ..., "-i", "--irma-directory", help="IRMA output directory to parse"
    ),
    jsonout: Path = typer.Option(..., "-j", "--json", help="json output path and name"),
    samplename: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name - included in the output json"
    ),
):
    """
    Parse IRMA output - return json
    """

    def get_irma_files(irmadir: Path, sample_name: str) -> dict:
        """
        Get IRMA output files of interest
        returns dict with absolute path

        *.coverage.a2m.txt
        *.vcf
        *.bam
        READ_COUNTS.txt
        """

        d = defaultdict()
        p = irmadir.resolve()
        # irma/sampleName/
        d["vcf"] = [str(x) for x in p.glob("*.vcf")]
        d["bam"] = [str(x) for x in p.glob("*.bam")]
        d["a2m"] = [str(x) for x in p.glob("*.a2m")]
        d["fasta"] = {x.stem: str(x) for x in p.glob("*.fasta")}

        # irma/sampleName/
        d["readcounts"] = (
            [str(p / "tables" / "READ_COUNTS.txt")]
            if (p / "tables" / "READ_COUNTS.txt").exists()
            else []
        )
        d["coverage"] = [str(x) for x in (p / "tables").glob("*coverage.a2m.txt")]
        d["insertions"] = [str(x) for x in (p / "tables").glob("*insertions.txt")]
        d["deletions"] = [str(x) for x in (p / "tables").glob("*deletions.txt")]
        d["variants"] = [str(x) for x in (p / "tables").glob("*variants.txt")]
        d["alleles"] = [str(x) for x in (p / "tables").glob("*allAlleles.txt")]
        # other
        d["secondary"] = {x.stem: str(x) for x in (p / "secondary").glob("*.fa")}
        d["amended"] = {
            x.suffix.replace(".", ""): str(x)
            for x in (p / "amended_consensus").glob("*")
        }

        d["runinfo"] = (
            str(p / "logs" / "run_info.txt")
            if (p / "logs" / "run_info.txt").exists()
            else []
        )

        # status represents the success or fail  or an irma run
        # at minimum a bam, fasta, vcf output is required
        if d["bam"] and d["fasta"] and d["vcf"]:
            d["status"] = "pass"
        else:
            d["status"] = "fail"

        final = {sample_name: d}
        return final

    parsed = get_irma_files(irmadir=irmadir, sample_name=samplename)

    with open(jsonout, "w", encoding="utf-8") as file:
        json.dump(parsed, file, ensure_ascii=False, indent=4)


@app.command(no_args_is_help=True)
def maskfasta(
    a2mfasta: Path = typer.Option(
        ..., "-a", "--a2m-fasta", help="path to the *.a2m  (fasta)"
    ),
    a2mcov: Path = typer.Option(
        ..., "-c", "--a2m-cov", help="The *coverage.a2m.txt (tsv) file from IRMA output"
    ),
    sample_name: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name, used to output"
    ),
    outdir: Path = typer.Option(..., "-o", "--output-dir", help="Output directory"),
    samcov: Path = typer.Option(
        None,
        "--samtools-depth",
        help="The samtools depth output - currently not supported",
    ),
    chromCol: int = typer.Option(
        0,
        "--chrom-col",
        help="Zero base position of the file column containing chromosome information",
    ),
    depthCol: int = typer.Option(
        6,
        "--depth-col",
        help="Zero base position of the file column containing depth information",
    ),
    sep: str = typer.Option(
        "\t", "--seperator", help="Seperator for out file \t or , "
    ),
    minCov: int = typer.Option(
        20,
        "-m",
        "--min-cov",
        help="Minimum coverage value. Any position with a depth below this value will be output",
    ),
):
    """
    Mask low coverage region using `bedtools maskfasta`

    irma settings: chromCol = 0 depthCol = 6

    samtools settings: % samtools depth -aa sample.sorted.bam > input.txt
    chromCol = 0 depthCol = 2
    """

    # outputs
    bed_out = outdir / f"{sample_name}.bed"
    masked_fasta = outdir / f"{sample_name}.masked.fasta"

    # read in depth file
    try:
        dat = pd.read_csv(a2mcov, sep=sep)
    except OSError as e:
        raise typer.BadParameter(f"Error reading infile file while trying to create bed file\n {e}")

    chrom = dat.iloc[:, chromCol][0]

    # select only Deletions and Matches
    dat = dat[dat["Alignment_State"].isin(["D", "M"])]

    # find regions less than the mincov
    s = pd.Series(dat.iloc[:, depthCol] <= minCov)
    grp = s.eq(False).cumsum()
    arr = (
        grp.loc[s.eq(True)].groupby(grp).apply(lambda x: [x.index.min(), x.index.max()])
    )

    # write out
    with open(outdir / f"{sample_name}.bed", "w", encoding="utf-8") as f:
        for l in list(arr):
            f.write(f"{str(chrom)}\t{l[0]}\t{l[1]}\n")

    # run bedtools
    subprocess.run(
        [
            "bedtools",
            "maskfasta",
            "-fi",
            a2mfasta,
            "-bed",
            bed_out,
            "-fo",
            masked_fasta,
        ],
        check=True,
    )

    if not masked_fasta.exists():
        raise typer.BadParameter("Error - no output was detected. Check inputs.")


@app.command(no_args_is_help=True)
def plot(
    vcffile: Path = typer.Option(
        ..., "-v", "--vcf-file", help="Path to the vcf"
    ),
    a2mcov: Path = typer.Option(
        ..., "-c", "--a2m-cov", help="Path to *coverage.a2m.txt (tsv)"
    ),
    sample_name: str = typer.Option(
        None, "-s", "--sample-name", help="Sample name, used to output"
    ),
    extra: str = typer.Option(
        None, "--extra", help="Extra information to include in the title"
    ),
    outfile: str = typer.Option(..., "-o", "--output-dir", help="Output file"),
    encoding: str = typer.Option(
        "svg", "-e", "--encoding", help="File encoding: svg, pdf, png, html"
    ),
):
    """
    Plot depth histograph annotated with potential ambiguous bases
    """

    a2m_df = pd.read_csv(
        a2mcov,
        sep="\t",
        dtype={
            "Reference_Name":"str",
            "Position":"float",
            "Coverage Depth":"int",
            "Consensus":"str",
            "Deletions":"int",
            "Ambiguous":"int",
            "Consensus_Count":"int",
            "Consensus_Average_Quality":"float",
            "HMM_Position":"int",
            "Alignment_state":"str"
            }
        )
    a2m_df = a2m_df[
        [
            'HMM_Position',
            'Position',
            'Coverage Depth',
            'Consensus',
            'Consensus_Count',
            'Alignment_State'
            ]
        ]
    vcf = pd.read_csv(
        vcffile,
        sep="\t",
        skiprows=21,
        names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    )

    vcf = vcf[["POS", "REF", "ALT", "QUAL", "FILTER", "INFO"]].rename(columns={
        "POS" : "Position",
        "REF" : "Reference",
        "ALT" : "Alternative",
        "QUAL" : "Quality",
        "FILTER" : "Filter",
        "INFO" : "Info"
    })
    vcf = pd.merge(vcf, a2m_df, how="left", on=['Position'])

    canvas = toyplot.Canvas(width=600, height=400)
    axes = canvas.cartesian(
        label=f"{sample_name} {extra or ''}",
        xlabel="HMM Position",
        ylabel="Depth"
    )

    covmax = max(a2m_df['Coverage Depth'])
    # toyplot
    mark = axes.fill(a2m_df['HMM_Position'], a2m_df['Coverage Depth'])
    axes.x.ticks.locator = toyplot.locator.Explicit(np.arange(0, 16000, 1000).tolist()) # type: ignore
    axes.x.ticks.show = True
    axes.y.ticks.show = True
    label_style = {"text-anchor":"start", "-toyplot-anchor-shift":"2.5px"}

    # add vertical lines with labels
    line = axes.vlines(vcf['HMM_Position'])
    for _, row in vcf.iterrows():
        axes.text(
            row['HMM_Position'],
            covmax + 10,
            row['Alternative'],
            style=label_style,
            color="red"
            )

    if encoding == "pdf":
        pdf.render(canvas, outfile)
    if encoding == "png":
        png.render(canvas, outfile)
    if encoding == "svg":
        svg.render(canvas, outfile)
    if encoding == "html":
        html.render(canvas, outfile)

@app.command(no_args_is_help=True,
             epilog="[yellow]Positions must be relative to the IRMA reference.[/yellow]")
def mutations(
    positions: str = typer.Option(
        None,
        "-p",
        "--positions",
        help='A list of positions to extract (comma separated and quoted). eg -p "31,34,140"',
    ),
    positions_file: Path = typer.Option(
        None,
        "-f",
        "--position-file",
        help="File containing positions to extract, one per line. ",
    ),
    alleles: Path = typer.Option(
        ..., "-a", "--alleles-file", help="IRMA *allAlleles.txt output"
    ),
    a2mcov: Path = typer.Option(
       ..., "-c", "--a2m-cov", help="The *coverage.a2m.txt (tsv) file from IRMA output"
    ),
    outfile: Path = typer.Option(..., "-o", "--outfile", help="Outfile (tsv)"),
):
    """
    !Not Implemented! Get mutations of interest from irma alleles.txt and vcf file.
    Positions
    """
    if positions and positions_file:
        raise typer.BadParameter(
            "--positions and --position-file are mutally exclusive. Use one, not both"
        )

    positions = [int(x) for x in positions.split(",")] # type: ignore

    try:
        allele_df = pd.read_csv(
            alleles,
            sep="\t",
            dtype={
                "Reference_Name":"str",
                "Position":"float",
                "Allele":"str",
                "Count":"int",
                "Total":"int",
                "Frequency":"float",
                "Average_Quality":"float",
                "ConfidenceNotMacErr":"float",
                "PairedUB":"float",
                "QualityUB":"float",
                "Allele_Type":"str",
                "HMM_Position":"int",
                "Alignment_state":"str"
                }
            )
        allele_df = allele_df.drop(["Reference_Name", "Position"], axis=1)
        a2m_df = pd.read_csv(
            a2mcov,
            sep="\t",
            dtype={
                "Reference_Name":"str",
                "Position":"float",
                "Coverage Depth":"int",
                "Consensus":"str",
                "Deletions":"int",
                "Ambiguous":"int",
                "Consensus_Count":"int",
                "Consensus_Average_Quality":"float",
                "HMM_Position":"int",
                "Alignment_state":"str"
                }
            )
    except OSError as e:
        raise typer.BadParameter(f"Error reading input files \n {e}")

    # HMM positions are used here as they're consist

    df = pd.merge(a2m_df, allele_df, how="left" ,on=["HMM_Position", "Alignment_State"])

    df = df[df['HMM_Position'].isin(positions)]
    df.to_csv(outfile, index=False, sep="\t")

    # if df.empty:
    #     raise typer.BadParameter("Empty Dataframe - Positions specified are not in Alleles file.\nThis occurs if the assembly is missing those positions")
    # print(df)

if __name__ == "__main__":
    app()