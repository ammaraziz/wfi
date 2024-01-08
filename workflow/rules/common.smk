'''
From the amazing https://github.com/gbouras13/hybracter
'''

import sys
import os
import pandas as pd
from  pathlib import Path
from collections import defaultdict

def get_input_r1(wildcards):
    return dictReads[wildcards.sample]["R1"]

def get_input_r2(wildcards):
    return dictReads[wildcards.sample]["R2"]

def samplesFromCsvShort(csvFile):
    """
    Read samples and files from a CSV Hybrid
    5 cols
    1 = sample
    2 = Long read
    3 = MinChromLength
    4 = R1 Short
    5 = R2 Short
    """
    outDict = {}
    with open(csvFile, "r") as csv:
        for line in csv:
            l = line.strip().split(",")
            if len(l) == 5:
                outDict[l[0]] = {}
                if (
                    os.path.isfile(l[1])
                    and l[2].isnumeric()
                    and os.path.isfile(l[3])
                    and os.path.isfile(l[4])
                ):
                    outDict[l[0]]["LR"] = l[1]
                    outDict[l[0]]["MinChromLength"] = l[2]
                    outDict[l[0]]["R1"] = l[3]
                    outDict[l[0]]["R2"] = l[4]
                else:
                    sys.stderr.write(
                        "\n"
                        f"    FATAL: Error parsing {csvFile}. One of \n"
                        f"    {l[1]} or \n"
                        f"    {l[3]} or \n"
                        f"    {l[4]} \n"
                        f"    does not exist or  {l[2]} is not an integer. \n"
                        "    Check formatting, and that \n"
                        "    file names and file paths are correct.\n"
                        "\n"
                    )
                    sys.exit(1)
            else:
                sys.stderr.write(
                    "\n"
                    f"    FATAL: Error parsing {csvFile}. Line {l} \n"
                    f"    does not have 5 columns. \n"
                    f"    Please check the formatting of {csvFile}. \n"
                )
                sys.exit(1)
    return outDict


def parseSamples(csvfile, long_flag):
    if os.path.isfile(csvfile) and long_flag is True:
        sampleDict = samplesFromCsvLong(csvfile)
    elif os.path.isfile(csvfile) and long_flag is False:
        sampleDict = samplesFromCsvShort(csvfile)
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

def getIRMAFiles(irma_output: str) -> dict[str]:
    '''
    Get IRMA output files of interest
    returns dict with absolute path

    *coverage.a2m.txt
    *.vcf
    *.bam
    READ_COUNTS.txt
    '''
    d = defaultdict(str)
    p = Path(irma_output)

    d['vcf'] = [str(x) for x in p.glob("*.vcf")]
    d['coverage'] = [str(x) for x in (p/"tables").glob("*coverage.a2m.txt")]
    d['bam'] = [str(x) for x in p.glob("*.bam")]
    d['counts'] = str(p/"tables"/"READ_COUNTS.txt") if (p/"tables"/"READ_COUNTS.txt").exists() else []

    return(d)

def make_bed_for_masking(a2m: str, min_cov:int=20) -> None:
    '''
    Given *coverage.a2m.txt input, create bed file with low coverage
    '''
    try:
        dat = pd.read_csv(a2m, sep='\t')
    except Exception as e:
        print(f"Error reading a2m file while trying to create bed file\n {e}")
    # find positions which are below cutoff
    dat = dat[dat['Alignment_State'].isin (['D', 'M'])]

    b = list(dat['Coverage Depth'] < 20)
    s = pd.Series(b)
    grp = s.eq(False).cumsum()
    arr = grp.loc[s.eq(True)] \
            .groupby(grp) \
            .apply(lambda x: [x.index.min(), x.index.max()])
    print(arr)
    return(arr)

#make_bed_for_masking("/Users/aaziz/repos/wfi/tests/outputirma/RSV321_S1/tables/rsv_a2-coverage.a2m.txt")