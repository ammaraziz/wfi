'''
Common.smk
'''

import sys
import os
from  pathlib import Path
from collections import defaultdict

import pandas as pd


def samplesFromCsv(csvFile):
    """
    Read samples and files from a CSV
    3 cols:
        1 = sampleid
        2 = runtype
        2 = R1 Short
        3 = R2 Short
    From the amazing https://github.com/gbouras13/hybracter
    """
    outDict = {}
    with open(csvFile, "r", encoding="utf-8") as csv:
        for line in csv:
            l = line.strip().split(",")
            if l[0].startswith("#"):
                continue
            if len(l) == 4:
                outDict[l[0]] = {}
                if (
                    #type(l[1]) is str
                    isinstance(l[1], str)
                    and os.path.isfile(l[2])
                    and os.path.isfile(l[3])
                ):
                    outDict[l[0]]["runtype"] = l[1]
                    outDict[l[0]]["R1"] = l[2]
                    outDict[l[0]]["R2"] = l[3]
                else:
                    sys.stderr.write(
                        "\n"
                        f"    Error parsing {csvFile}.\n"
                        f"    {l[1]} must be paired-end or single-end or long or \n"
                        f"    {l[2]} must exist or \n"
                        f"    {l[3]} must exist or \n"
                        "    Check formatting, and that \n"
                        "    file names and file paths are correct.\n"
                        "\n"
                    )
                    sys.exit(1)
            else:
                sys.stderr.write(
                    "\n"
                    f"    FATAL: Error parsing {csvFile}. Line {l} \n"
                    f"    does not have 4 columns. \n"
                    f"    Please check the formatting of {csvFile}. \n"
                )
                sys.exit(1)
    return outDict


def parseSamples(csvfile):
    """
    From the amazing https://github.com/gbouras13/hybracter
    """
    if os.path.isfile(csvfile):
        sampleDict = samplesFromCsv(csvfile)
    else:
        sys.stderr.write(
            "\n"
            f"    FATAL: something is wrong. Likely {csvfile} is neither a file nor directory.\n"
            "\n"
        )
        sys.exit(1)

    # checks for dupes
    SAMPLES = list(sampleDict.keys())

    # Check for duplicates
    has_duplicates = len(SAMPLES) != len(set(SAMPLES))

    # error out if dupes
    if has_duplicates is True:
        sys.stderr.write(
            f"Duplicates found in the SAMPLES list in column 1 of {csvfile}.\n"
            f"Please check {csvfile} and give each sample a unique name!"
        )
        sys.exit(1)

    return sampleDict

# get inputs
def get_input_r1(wildcards):
    return dictReads[wildcards.sample]["R1"]

def get_input_r2(wildcards):
    return dictReads[wildcards.sample]["R2"]

def get_input_lr_fastqs(wildcards):
    return dictReads[wildcards.sample]["LR"]

def get_irma_files(irma_output: str) -> dict:
    '''
    Get IRMA output files of interest
    returns dict with absolute path

    *.coverage.a2m.txt
    *.vcf
    *.bam
    READ_COUNTS.txt
    '''
    d = defaultdict()
    p = Path(irma_output)

    d['vcf'] = [str(x) for x in p.glob("*.vcf")]
    d['coverage'] = [str(x) for x in (p/"tables").glob("*coverage.a2m.txt")]
    d['bam'] = [str(x) for x in p.glob("*.bam")]
    d['counts'] = str(p/"tables"/"READ_COUNTS.txt") if (p/"tables"/"READ_COUNTS.txt").exists() else []

    return d

def make_bed_for_masking(infile: str,
                         chromCol: int,
                         depthCol: int,
                         outfile: str,
                         sep: str="\t",
                         min_cov:int=20) -> None:
    '''
    Create bed file from tsv file containing depths.
    
    chromCol:   pos of chromosome column.
    depthCol:   pos of column with coverage values.
    all are 0-based positions.

    irma settings:
    *-coverage.a2m.txt
    chromCol = 0
    depthCol = 6

    samtools settings:
    samtools depth -aa sample.sorted.bam > depth.txt
    chromCol = 0
    depthCol = 2
    '''

    try:
        dat = pd.read_csv(infile, sep=sep)
    except OSError as e:
        print(f"Error reading infile file while trying to create bed file\n {e}")
        sys.exit(1)
    chrom = dat.iloc[:, chromCol][0]

    # irma specific 
    dat = dat[dat['Alignment_State'].isin (['D', 'M'])]

    s = pd.Series(dat.iloc[:, depthCol] < min_cov)
    grp = s.eq(False).cumsum()
    arr = grp.loc[s.eq(True)] \
            .groupby(grp) \
            .apply(lambda x: [x.index.min(), x.index.max()])

    with open(outfile, "w", encoding="utf-8") as f:
        for l in list(arr):
            f.write(f"{str(chrom)}\t{l[0]}\t{l[1]}\n")

    return arr

# make_bed_for_masking("/Users/aaziz/repos/wfi/tests/output/irma/RSV321_S1_L001/tables/rsv_a2-coverage.a2m.txt", 
# chromCol=0,
# depthCol=6,
# outfile="/Users/aaziz/repos/wfi/tests/output/irma/mask.bed")

def get_irma_module_config(org: str, technology: str, secondary: bool = False) -> str:
    """
    Return correct irma module for the organism/technology
    """
    org = org.lower()
    technology = technology.lower()

    if org == "rsv" and bool is False:
        raise KeyError("IRMA RSV module does not support secondary assembly")

    modules = {
        "flu" : {
            "illumina" : "FLU",
            "nanopore" : "FLU-minion",
            "pgm" : "FLU-pgm"
            },
        "rsv" : {
            "illumina" : "RSV",
            "nanopore" : "RSV-minion",
        },
    }
    if secondary:
        m = modules[org][technology] + "-secondary"
    else:
        m = modules[org][technology]
    return m
