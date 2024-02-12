"""
IRMA Helpers - parse output directory
"""
import json
import subprocess
from math import ceil
from pathlib import Path
from collections import defaultdict
import random

import prepare

from Bio import SeqIO

import typer
from rich import print as pprint
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
    irmadir: Path = typer.Option(..., "-i", "--irma-dir", help="IRMA output directory to parse"),
    jsonout: Path = typer.Option(..., "-j", "--json", help="json output path and name"),
    samplename: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name - included in the output json"
    ),
    organism: str = typer.Option(
        ..., "--organism", help="IRMA Module / Organism [RSV, FLU, FLU_AD, EBOLA, COV]"
    )
):
    """
    Parse IRMA output - return json.
    """

    def get_irma_files(irmadir: Path, samplename: str) -> dict:
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
            else [])[0]

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

        # info
        d["runinfo"] = (
            str(p / "logs" / "run_info.txt")
            if (p / "logs" / "run_info.txt").exists()
            else []
        )

        d["masked"] = {}
        d["depth_plots"] = {}
        d["statistics"] = {}
        # status represents the success or fail  or an irma run
        # at minimum a bam, fasta, vcf output is required
        if d["bam"] and d["fasta"] and d["vcf"]:
            d["status"] = "pass"
        else:
            d["status"] = "fail"

        return {samplename: d}

    # genotype
    def genotype(a2m: dict, organism: str) -> tuple[str, str]:
        '''
        Extract subtype from IRMA named output of a2m
        return two values:
            - species (RSV, FLUA, FLUB... etc, EBOLA)
            - subtype if applicable, eg for FLU HxNx
        '''

        # for RSV items is ['org', 'subtype']
        # for FLU HA/NA items is ['org', 'gene', 'subtype']
        # for FLU other genes ['org', 'gene']
        # for all others ['org']
        if organism in ["FLU", "FLU-AD"]:
            FLU = None
            HA = None
            NA = None

            for name in a2m.keys():
                parts = sorted(name.split("_"))

                if parts[1] in "HA": # envlope genes
                    HA = f"H{parts[2]}"
                    FLU = f"FLU{parts[0]}"
                if parts[1] in "NA": # envlope genes
                    NA = f"N{parts[2]}"
                    FLU = f"FLU{parts[0]}"
                if parts[1] in ["PB2","PB1","PA","NP","MP","NS"]:
                    FLU = f"FLU{parts[0]}"

            if HA is None:
                HA = "Hx"
            if NA is None:
                NA = "Nx"
            if FLU is None: # catch all, shoudn't happen
                FLU = "FLUx"

            return FLU, f"{HA}{NA}"
        for k in a2m.keys():
            parts = k.split("_")
            if organism == "RSV":
                ORG = parts[0].upper()
                TYPE = parts[1]
                return ORG, TYPE

        # all other organisms that do not have subtypes
        parts = list(a2m.keys())[0].split("_")
        ORG = parts[0].upper()
        TYPE = "NA"
        return ORG, TYPE

    parsed = get_irma_files(irmadir=irmadir, samplename=samplename)
    org, gentype = genotype(parsed[samplename]['a2m'], organism=organism)
    parsed[samplename]['organism'] = org
    parsed[samplename]['type'] = gentype

    with open(jsonout, "w", encoding="utf-8") as file:
        json.dump(parsed, file, ensure_ascii=False, indent=4)

@app.command(no_args_is_help=True, epilog="[yellow]Note: .plot.json is used[/yellow]")
def cjson(
    pipeline_outdir: Path = typer.Option(
        ..., "--pipeline-outdir", help="wfi output directory containing /irma/[sampleName].json"
    ),
    pattern: str = typer.Option(
        ".plot.json", "--pattern", help="pattern to match json output to merge. Do not include regex pattern."
    ),
    cjson_out: Path = typer.Option(
        ..., "--json-out", help="Output combined Json file"
    ),
    ):
    '''
    Combine the [sampleName].json of all samples in a run.
    '''
    jfiles = [x for x in pipeline_outdir.glob(f"**/*{pattern}")]

    cjson = {}
    for j in jfiles:
        with open(j, "r") as handle:
            temp_json = json.loads(handle.read())
            cjson = cjson | temp_json

    with open(cjson_out, "w") as handle:
        json.dump(cjson, handle, ensure_ascii=False, indent=4)
        pprint(f"[green]Combined json written to: {cjson_out.absolute()}[/green]")

@app.command(
        no_args_is_help=True,
        epilog="[yellow]Note: json file is ammended in place to include new giles generated.[/yellow]")
def maskfasta(
    json_file: Path = typer.Option(
        ..., "-j", "--json", help="irmakit.py parse output"
    ),
    sample_name: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name, used to output"
    ),
    json_out: Path = typer.Option(
        ..., "--json-out", help="Amended json file output"
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
    Mask low coverage region using `bedtools maskfasta`.

    irma settings: chromCol = 0 depthCol = 6
    """

    outdir.mkdir(exist_ok=True, parents=True)

    # read in depth file
    try:
        with open(json_file, "r", encoding="utf-8") as handle:
            irma_data = json.loads(handle.read())
            a2m_cov = irma_data[sample_name]['coverage']
            a2m_fasta = irma_data[sample_name]['a2m']
    except Exception as e:
        raise typer.BadParameter(f"Could not read json input.\n{e}")

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

        # ammend json
        irma_data[sample_name]['masked'].update({fname : {"bed" : str(bed_out), "masked_fasta" : str(masked_fasta)}})

        with open(json_out, "w", encoding="utf-8") as handle:
            json.dump(irma_data, handle, ensure_ascii=False, indent=4)

        if not masked_fasta.exists():
            raise typer.BadParameter("Error - no output was detected. Check inputs.")

    pprint(f"[green]Masking complete: {str(masked_fasta.resolve().parent)}[/green]") # type: ignore

@app.command(no_args_is_help=True)
def plot(
    json_file: Path = typer.Option(
        ..., "-j", "--json", help="irmakit.py parse output"
    ),
    sample_name: str = typer.Option(
        ..., "-s", "--sample-name", help="Sample name, used to output"
    ),
    outfile: Path = typer.Option(..., "-o", "--output-file", help="Output file"),
    json_out: Path = typer.Option(
        ..., "--json-out", help="Amended json file output"
    ),
    extra: str = typer.Option(
        None, "--extra", help="Extra information to include in the title"
    ),
    encoding: str = typer.Option(
        "html", "-e", "--encoding", help="File encoding: svg, pdf, png, html"
    ),
):
    """
    Plot depth histograph annotated with potential ambiguous bases.
    """

    if encoding not in ['svg', 'pdf', 'png', 'html']:
        raise typer.BadParameter(f"Encoding {encoding} not supported. Possible values: svg, pdf, png, html")

    with open(json_file, "r", encoding="utf-8") as handle:
        irma_data = json.loads(handle.read())
        a2m_cov = irma_data[sample_name]['coverage']
        vcf_file = irma_data[sample_name]['vcf']
        masked = irma_data[sample_name]['masked']

    datasets = defaultdict(dict)
    for fname in a2m_cov.keys():

        # it's assumed for each assembled gene/organism an coverage-a2m file is created
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
        # vcf files
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

        bed = pd.read_csv(masked[fname]['bed'], sep="\t",
                           names=["chrom", "start", "end"])

        ### compute data for plotting
        misc_size = -10
        genome_size = a2m_df.shape[0]
        # gap column plotting
        a2m_df['Gaps'] = np.where(a2m_df['Alignment_State'] != 'D', 0,
                                np.where(a2m_df['Alignment_State'] == 'D', misc_size, 0))
        # masked positions
        np_masked_region = np.zeros(genome_size)
        for index, row in bed.iterrows():
            np_masked_region[row['start'] : row['end']] = misc_size
        a2m_df['Masked Regions'] = np_masked_region

        # output data for plotting
        datasets[fname] = {
            "vcf" : vcf,
            "a2m_df" : a2m_df,
            "covmax" : max(a2m_df['Coverage Depth']),
            "genomesize_upper" : ceil(genome_size / 1000) * 1000, # round UP to the nearest 1000
            "interval" : round(a2m_df.shape[0], -3)/10,
            }

    # plot data
    canvas = toyplot.Canvas(width=1000, height=800)
    for index, (fname, data) in enumerate(datasets.items()):

        axes = canvas.cartesian(
            label=f"{sample_name} - {fname} {extra or ''}",
            xlabel="HMM Position",
            ylabel="Depth",
            grid=(len(datasets), 2, index)
            )
        # toyplot
        depth = axes.fill(data['a2m_df']['HMM_Position'], data['a2m_df']['Coverage Depth'], title="Depth")
        masked = axes.fill(data['a2m_df']['HMM_Position'], data['a2m_df']['Masked Regions'], color = "orange", title="Masked Regions")
        gaps = axes.fill(data['a2m_df']['HMM_Position'], data['a2m_df']['Gaps'], color = "red", title = "Regions without Coverage")
        # ticks
        axes.x.ticks.locator = toyplot.locator.Explicit( # type: ignore
            np.arange(0, data['genomesize_upper'], data['interval']).tolist()
            )
        axes.x.ticks.show = True
        axes.y.ticks.show = True
        label_style = {"text-anchor":"start", "-toyplot-anchor-shift":"2px", "font-size" : 10}

        # display mask cutoff
        axes.hlines(20, style={"stroke":"blue", "stroke-dasharray":"2, 2"})
        axes.text(genome_size/2, 25, "Cutoff Masking", color = "black", style={"font-size" : 10}) # type: ignore

        # add vertical lines with labels
        axes.vlines(data['vcf']['HMM_Position'], title="Ambig Bases")
        for index, row in data['vcf'].iterrows():
            axes.text(
                row['HMM_Position'],
                data['covmax'] + random.randint(1,10),
                row['Alternative'],
                style=label_style,
                color="red",
                title="Ambiguous Bases"
                )
        # legend
        canvas.legend(
            [("Gaps", gaps),
             ("Masked Regions", masked)],
            corner=("top-right", 50, 100, 35),
        )

    eval(encoding).render(canvas, str(outfile.resolve()))

    # ammend json
    irma_data[sample_name]['depth_plots'].update({encoding : str(outfile.resolve())})

    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(irma_data, handle, ensure_ascii=False, indent=4)

        pprint(f"[green]Plotting complete - {len(datasets)} graphs in {outfile.resolve()} [/green]")


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
    json_out: Path = typer.Option(None, "--json-out", help="Amended json file output"),
):
    """
    !Not Implemented! Get mutations of interest from irma alleles.txt and vcf file.
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


@app.command(no_args_is_help=True)
def stats(
    json_file: Path = typer.Option(
        ..., "-j", "--json", help="irmakit.py parse output"
    ),
    json_out: Path = typer.Option(
        ..., "--json-out", help="Amended json file output"
    ),
    ):
    '''
    Per run (irma) statistics. Json input is modified inplace. 
    '''
    with open(json_file, "r", encoding="utf-8") as handle:
        irma_data = json.loads(handle.read())
        sample_name = next(iter(irma_data))
        read_counts = pd.read_csv(irma_data[sample_name]['readcounts'], sep="\t", index_col="Record")

    mean_depth = {}
    num_gaps = {}

    for key in irma_data[sample_name]['a2m'].keys():
        dat = pd.read_csv(irma_data[sample_name]['coverage'][key], sep="\t")
        mean_depth[key] = round(dat.loc[:, "Coverage Depth"].mean(), 2)
        num_gaps[key] = int(dat.loc[:, 'Alignment_State'].value_counts()['D'])

    irma_data[sample_name]['statistics'].update(
        {
            "read_counts" : {
                "initial" : int(read_counts.loc['1-initial', 'Reads']), # type: ignore
                "passqc" : int(read_counts.loc['2-passQC', 'Reads']), # type: ignore
                "match" : int(read_counts.loc['3-match', 'Reads']), # type: ignore
                "altmatch" : int(read_counts.loc['3-altmatch', 'Reads']), # type: ignore
                "patterm_match" : read_counts.iloc[9:, 1].to_dict(),
                },
            "mean_depth" : mean_depth,
            "num_gaps" : num_gaps,
        }
    )
    with open(json_out, "w", encoding="utf-8") as handle:
        json.dump(irma_data, handle, ensure_ascii=False, indent=4)

@app.command(no_args_is_help=True,
        epilog="[yellow]Note: Combine is only for segemented organisms such as FLU")
def crs(
    cjson_file: Path = typer.Option(..., "--cjson", help="Combined json file."),
    outdir: Path = typer.Option(..., "--outdir", help="output directory"),
    ):
    """
    Combine -> Rename -> Sort IRMA output
    """

    segements = {
        # Flu A/B
        "PB2": "1",
        "PB1": "2",
        "PA": "3",
        "HA": "4",
        "NP": "5",
        "NA": "6",
        "MP": "7",
        "NS": "8",
        # Flu C
        "HE": "4",
        "P3": "3",
        # RSV
        "A1": "A1",
        "A2": "A2",
        "A3": "A2",
        "B": "B",
    }

    with open(cjson_file, "r", encoding="utf-8") as handle:
        cjson = json.loads(handle.read())

    # read in all sequences into a dict of lists for each IRMA run
    sequences = defaultdict(list)
    for sample,contents in cjson.items():
        sequences[sample] = []
        for _,p in contents['a2m'].items():
            for record in SeqIO.parse(p, "fasta"):
                sequences[sample].append(record)

    # for each sample, rename + combine sequences
    for sample,contents in sequences.items():
        for record in contents:
            # rename
            seg = segements[record.id.upper().split("_")[1]]
            record.id = f"{sample}.{seg}" #segements[id.split("_")[1]
            record.description = ""
        # combine
        with open(outdir / f"{sample}.fasta", "w", encoding="utf-8") as handle:
            SeqIO.write(contents, handle, "fasta")

@app.callback()
def callback():
    """
    Typer app, including Click subapp
    """

typer_click_object = typer.main.get_command(app)

typer_click_object.add_command(prepare.prepare) # type: ignore

if __name__ == "__main__":
    typer_click_object()
