#!/usr/bin/env python

"""
Process and subtype IRMA output
subtyper.py -i [inputDir] -o [outputdir] -r [species] --subset [all, subset rsv]
inputDir is expected to be the /irma/ folder
"""

import shutil
import os
import re
import argparse
from pathlib import Path
from Bio import SeqIO


def getSegment(header: str) -> int:
    genes = {
        "PB2": 1,
        "PB1": 2,
        "PA": 3,
        "HA": 4,
        "NP": 5,
        "NA": 6,
        "MP": 7,
        "NS": 8,
        # Flu C
        "HE": 4,
        "P3": 3,
        # RSV
        "A1": 1,
        "A2": 1,
        "B": 1,
    }
    return(genes[header.upper()])


def getSamples(inputFolder: str) -> list:
    """
    Find all IRMA output from a [inputFolder] returning complete path

    input
        inputFolder  :   str: directory path
    return
        list of file path objects
    """

    abs_path = [
        os.path.abspath(os.path.join(inputFolder, p))
        for p in os.listdir(inputFolder)
        if p not in ["rename", "consensus", "bySubtype"] and not p.startswith("NTC_")
    ]
    output = list(filter(os.path.isdir, abs_path))
    return output


def combine_fasta(directory, subset, outname="all.fasta"):
    """
    combine all IRMA output or specific subset in a directory

    subset is segments to keep:
    all = all segments
    subsetA = ha, na, mp
    subsetB = ha, na

    output = all.fasta
    return - Path of all.fasta
    """

    subsetDict = {
        "all": ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"],
        "subset": ["HA", "NA", "MP"],
        "rsv": ["a1", "a2", "b"],
    }

    fastaToCombine = [
        a
        for a in list(Path(str(directory)).glob("*.fasta"))
        if os.path.basename(a) not in ["contigs.fasta", "all.fasta"]
        if a.stem.split("_")[1] in subsetDict[subset]
    ]
    outpath = Path(directory, outname)

    with open(outpath, "w") as outfile:
        if not fastaToCombine:
            outfile.write(
                f"No fasta was found. IRMA did not assembly any viruses. Perform QC on: {directory}"
            )
            return outpath
        else:
            for filename in fastaToCombine:
                if filename == outpath:
                    continue
                else:
                    with open(filename, "r") as readfile:
                        outfile.write(readfile.read())
    return outpath


def get_subtype(sampleDirectory, org: str):
    """
    get subtype from output of IRMA file by parsing fasta file names
    """
    if org == "FLU":
        regex_spec = re.compile(r"([A|B])_\w+\.fasta")
        regex_ha = re.compile(r"[A|B]_HA_(H\d+)\.fasta")
        regex_na = re.compile(r"[A|B]_NA_(N\d+)\.fasta")
        for d in os.listdir(sampleDirectory):
            if regex_ha.search(d):
                ha = list(regex_ha.search(d).groups(1))
            if regex_na.search(d):
                na = list(regex_na.search(d).groups(1))
            if regex_spec.search(d):
                spec = list(regex_spec.search(d).groups(1))

        if "na" not in vars():
            na = [""]
        if "ha" not in vars():
            ha = [""]
        if "spec" not in vars():
            spec = [""]
        return spec + ha + na

    if org == "RSV":
        regex_spec = re.compile(r"rsv_(.+)\.fasta")
        for d in os.listdir(sampleDirectory):
            if regex_spec.search(d):
                spec = list(regex_spec.search(d).groups(1))
        if "spec" not in vars():
            regex_spec = [""]
        return [x.upper() for x in spec] + [""] + [""]


def makeSubtypeFolderStruct(outputDir: str, subtypeList: list) -> Path:
    """
    Given a list: [virus, ha, na] return a Path object that reflects output
    For FLU = ['A', 'H1', 'N2] returns Path(outputDir/A/H1N1)
    For RSV = ['A2'] = Path(outputDir/A2/)
    """
    return Path(outputDir).joinpath(subtypeList[0], subtypeList[1] + subtypeList[2])


    def renameFastaRecord(fastaRecordId: str, sampleName: str) -> str:
        """
        rename fasta file from IRMA output (eg A_HA_H3.fasta) to [Isolate].[segement_number]
        """
        segment = getSegment(fastaRecordId.split("_")[1])
        fastaRecordId = sampleName + "." + str(segment)

        return fastaRecordId


def register_arguments():
    parser = argparse.ArgumentParser(
        description="Combine, rename, subtype IRMA outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="subtyper.py -i [inputDir] -r [species] --subet []",
    )

    parser.add_argument(
        "-i", "--input", required=True, help="Input directory eg '/irma/"
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument(
        "-r", "--species", required=True, help="Virus species: flu, rsv"
    )
    parser.add_argument(
        "--subset",
        required=False,
        help="Subset output to specific gene segments: all, subsetA, subsetB, rsv",
        default="all",
    )
    return parser.parse_args()


def run() -> None:
    args = register_arguments()
    args.species = args.species.upper()

    out_renamed = Path(args.output).joinpath("consensus")
    out_subtyped = Path(args.output).joinpath("bySubtype")
    out_results = Path(args.output).joinpath("subtypes.tsv")
    os.makedirs(out_renamed, exist_ok=True)
    os.makedirs(out_subtyped, exist_ok=True)

    irma_output = getSamples(args.input)

    for folder in irma_output:
        sampleName = os.path.basename(folder).split("_")[0]

        # combine fasta
        allFastaLoc = combine_fasta(folder, args.subset)

        # rename headers and save in output folder
        allFastaRenamed = out_renamed.joinpath(sampleName + ".fasta")

        # catch no result/ empty file
        tmp = SeqIO.parse(allFastaLoc, "fasta")
        if not [len(record) for record in tmp]:
            continue

        # rename records
        with open(allFastaLoc, "r") as original, open(allFastaRenamed, "w") as renamed:
            for record in SeqIO.parse(original, "fasta"):
                record.id = renameFastaRecord(record.id, sampleName)
                record.description = record.id
                SeqIO.write(record, renamed, "fasta")

        # subtyping
        subtype = get_subtype(folder, args.species)
        with open(out_results, "a+") as results:
            results.write(str(Path(folder).name) + "\t" + "".join(subtype) + "\n")

        # create subtype folder structure
        outSubtypeLoc = makeSubtypeFolderStruct(out_subtyped, subtype)
        os.makedirs(outSubtypeLoc, exist_ok=True)

        # copy fasta file
        shutil.copy(allFastaRenamed, outSubtypeLoc.joinpath(allFastaRenamed.name))


if __name__ == "__main__":
    run()
