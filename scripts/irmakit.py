"""
IRMA Helpers - parse output directory
"""
import json
import subprocess
from math import ceil
from pathlib import Path
from collections import defaultdict

import typer
from rich import print
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
        # irma/sampleName/ - multiple files are possible
        d["vcf"] = {x.stem: str(x) for x in p.glob("*.vcf")}
        d["bam"] = {x.stem: str(x) for x in p.glob("*.bam")}
        d["a2m"] = {x.stem: str(x) for x in p.glob("*.a2m")}
        d["fasta"] = {x.stem: str(x) for x in p.glob("*.fasta")}

        # irma/sampleName/ - single file only
        d["readcounts"] = (
            [str(p / "tables" / "READ_COUNTS.txt")]
            if (p / "tables" / "READ_COUNTS.txt").exists()
            else []
        )[0]
        # multiple possible
        d["coverage"] = {x.stem.replace("-coverage.a2m", ""): str(x) for x in (p / "tables").glob("*coverage.a2m.txt")}
        # single instances only
        d["insertions"] = [str(x) for x in (p / "tables").glob("*insertions.txt")][0]
        d["deletions"] = [str(x) for x in (p / "tables").glob("*deletions.txt")][0]
        d["variants"] = [str(x) for x in (p / "tables").glob("*variants.txt")][0]
        d["alleles"] = [str(x) for x in (p / "tables").glob("*allAlleles.txt")][0]
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
    json_file: Path = typer.Option(
        ..., "-j", "--json", help="irmakit.py parse output"
    ),
    sample_name: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name, used to output"
    ),
    outdir: Path = typer.Option(..., "-o", "--output-dir", help="Output directory"),
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

    outdir.mkdir(exist_ok=True, parents=True)

    # read in depth file
    try:
        with open(json_file) as handle:
            irma_data = json.loads(handle.read())
            a2m_cov = irma_data[sample_name]['coverage']
            a2m_fasta = irma_data[sample_name]['a2m']
    except Exception as e:
        raise typer.BadParameter(f"Could not read json input. \n {e}")

    # ensure the same number of keys (files) are in both dicts
    if not len(a2m_cov.keys()) == len(a2m_fasta.keys()):
        raise typer.BadParameter(f"Number of files differ between coverage-a2m and fasta-a2m for {sample_name}")

    for fname in a2m_cov.keys():
        # outputs
        bed_out = outdir / f"{sample_name}_{fname}.bed"
        masked_fasta = outdir / f"{sample_name}_{fname}.masked.fasta"

        dat = pd.read_csv(a2m_cov[fname], sep="\t")
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
        with open(bed_out, "w", encoding="utf-8") as f:
            for l in list(arr):
                f.write(f"{str(chrom)}\t{l[0]}\t{l[1]}\n")

        # run bedtools
        subprocess.run(
            [
                "bedtools",
                "maskfasta",
                "-fi",
                a2m_fasta[fname],
                "-bed",
                bed_out,
                "-fo",
                masked_fasta,
            ],
            check=True,
        )

        if not masked_fasta.exists():
            raise typer.BadParameter("Error - no output was detected. Check inputs.")

    print(f"[green]Masking complete - Outputs: \n - {masked_fasta}\n - {bed_out}[/green]")


@app.command(no_args_is_help=True)
def plot(
    json_file: Path = typer.Option(
        ..., "-j", "--json", help="irmakit.py parse output"
    ),
    sample_name: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name, used to output"
    ),
    outfile: str = typer.Option(..., "-o", "--output-file", help="Output file"),
    extra: str = typer.Option(
        None, "--extra", help="Extra information to include in the title"
    ),
    encoding: str = typer.Option(
        "svg", "-e", "--encoding", help="File encoding: svg, pdf, png, html"
    ),
):
    """
    Plot depth histograph annotated with potential ambiguous bases
    """

    with open(json_file) as handle:
        irma_data = json.loads(handle.read())
        a2m_cov = irma_data[sample_name]['coverage']
        vcf_file = irma_data[sample_name]['vcf']

    datasets = defaultdict(dict)
    for fname in a2m_cov.keys():
        a2m_df = pd.read_csv(
            a2m_cov[fname],
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
            vcf_file[fname],
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
        datasets[fname] = {
            "vcf" : vcf,
            "a2m_df" : a2m_df,
            "covmax" : max(a2m_df['Coverage Depth']),
            "genomesize" : ceil(a2m_df.shape[0] / 1000) * 1000, # round UP to the nearest 1000
            "interval" : round(a2m_df.shape[0], -3)/10,
            }
    canvas = toyplot.Canvas(width=800, height=600)
    for index, (fname, data) in enumerate(datasets.items()):

        axes = canvas.cartesian(
            label=f"{sample_name} - {fname}",
            xlabel="HMM Position",
            ylabel="Depth",
            grid=(len(datasets), 2, index)
        )

        # toyplot
        axes.fill(data['a2m_df']['HMM_Position'], data['a2m_df']['Coverage Depth'])
        axes.x.ticks.locator = toyplot.locator.Explicit(np.arange(0, data['genomesize'], data['interval']).tolist()) # type: ignore
        axes.x.ticks.show = True
        axes.y.ticks.show = True
        label_style = {"text-anchor":"start", "-toyplot-anchor-shift":"2.5px"}

        # add vertical lines with labels
        line = axes.vlines(data['vcf']['HMM_Position'])
        for _, row in data['vcf'].iterrows():
            axes.text(
                row['HMM_Position'],
                data['covmax'] + 10,
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

    print(f"[green] Plotting complete - {len(datasets)} graphs in {outfile} [/green]")


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
