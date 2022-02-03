#!/usr/bin/env python

'''
Process and subtype IRMA output - flu only
'''

import shutil
import os
import re
import argparse
from pathlib import Path
from Bio import SeqIO

def getSegment(header):
    genes = {
        'PB2' : 1,
        'PB1' : 2,
        'PA'  : 3,
        'HA'  : 4,
        'NP'  : 5,
        'NA'  : 6,
        'MP'  : 7,
        'NS'  : 8,
        # Flu C
        'HE'  : 4,
        'P3'  : 3
        }
    try:
        seg = genes[header.upper()]
    except:
        raise ValueError("Unknown segment. Check input: " + str(header))
    return(seg)

def getSamples(inputFolder):
    '''
    Find all IRMA output from a [inputFolder] returning complete path

    input
        inputFolder  :   str: directory path
    return
        list of file path objects
    '''
    abs_path = [
        os.path.abspath(os.path.join(inputFolder, p)) for p 
        in os.listdir(inputFolder)
        if p not in ['rename', 'renamed', 'bySubtype'] and not p.startswith("NTC_")
        ]
    output = list(filter(os.path.isdir, abs_path))
    return(output)

def combineFasta(directory, subset, outname="all.fasta"):
    '''
    combine all IRMA output or specific subset in a directory 

    subset is segments to keep:
    all = all segments
    subsetA = ha, na, mp
    subsetB = ha, na

    output = all.fasta
    return - Path of all.fasta
    '''

    subsetDict = {
        'all' : ['PB2', 'PB1', 'PA', 'HA', ' NP', 'NA', 'MP', 'NS'],
        'subsetA' : ['HA', 'NA', 'MP'],
        'subsetB' : ['HA', 'NA']
    }
    segsToKeep = subsetDict[subset]


    fastaToCombine = [
        a for a 
        in list(Path(str(directory)).glob("*.fasta"))
        if os.path.basename(a) not in ['contigs.fasta', 'all.fasta']
        if a.stem.split('_')[1] in segsToKeep
    ]

    # check - fasta found?
    if not fastaToCombine:
        raise ValueError("Input directory did not contain any folders with Fasta files. Check input!")

    outpath = Path(directory, outname)

    with open(outpath, 'w') as outfile:
        if not fastaToCombine:
            outfile.write("No fasta was found. IRMA did not assembly any viruses. Perform QC.")
            return
        else:
            for filename in fastaToCombine:
                if filename == outpath:
                    continue
                else:
                    with open(filename, 'r') as readfile:
                        outfile.write(readfile.read())
    return(outpath)


def getSubtype(sampleDirectory):
    '''
    get subtype from output of IRMA file by parsing HA and NA fasta files
    '''

    # species
    regex_spec = re.compile(r'([A|B])_\w+\.fasta')

    regex_ha = re.compile(r'[A|B]_HA_(H\d+)\.fasta')
    regex_na = re.compile(r'[A|B]_NA_(N\d+)\.fasta')
    for d in os.listdir(sampleDirectory):
        if regex_ha.search(d):
            ha = list(regex_ha.search(d).groups(1))
        if regex_na.search(d):
            na = list(regex_na.search(d).groups(1))
        if regex_spec.search(d):
            spec = list(regex_spec.search(d).groups(1))
    
    if 'na' not in vars():
        na = ['']
    if 'ha' not in vars():
        ha = ['']
    if 'spec' not in vars():
        spec = ['']
    return(spec + ha + na)

def makeSubtypeFolderStruct(outputDir, subtypeList):
    '''
    Given a list: [virus, ha, na] return a Path object that reflects output
    For example ['A', 'H1', 'N2] returns Path(outputDir/A/H1N1) 
    '''
    return(Path(outputDir).joinpath(subtypeList[0], subtypeList[1]+subtypeList[2]))

def renameFastaRecord(fastaRecordId, sampleName):
    '''
    rename fasta file from IRMA output (eg A_HA_H3.fasta) to [Isolate].[segement_number]
    '''
    segment = getSegment(fastaRecordId.split("_")[1])
    fastaRecordId = sampleName + "." + str(segment)
    
    return(fastaRecordId)

def register_arguments():
    parser = argparse.ArgumentParser(
        description="Combine, rename, subtype IRMA outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog= "subtyper.py -i [inputDir] -o [outputDir] -r [species]"
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input directory eg '/assemblies/")
    # parser.add_argument("-o", "--output", required=True,
    #                     help="Output directory")
    parser.add_argument("-r", "--species", required=True,
                        help="Virus species: flu, rsv"),
    parser.add_argument("--subset", required=True,
                        help ="Subset output to specific gene segments: all, subsetA, subsetB, rsv",
                        default='all')
    return(parser.parse_args())

def run():
    args = register_arguments()

    out_renamed = Path(args.input).parent.joinpath("renamed")
    out_subtyped = Path(args.input).parent.joinpath('bySubtype')
    out_results = Path(out_subtyped / "subtypes.tsv")
    os.makedirs(out_renamed, exist_ok=True)
    os.makedirs(out_subtyped, exist_ok=True)
    
    irma_output = getSamples(args.input)
    
    for folder in irma_output:
        sampleName = os.path.basename(folder).split('_')[0]
        
        # combine fasta
        allFastaLoc = combineFasta(folder, args.subset)
        
        # rename headers and save in output folder
        allFastaRenamed = out_renamed.joinpath(sampleName + ".fasta")
        with open(allFastaLoc, 'r') as original, open(allFastaRenamed, 'w') as renamed:
            for record in SeqIO.parse(original, 'fasta'):
                record.id = renameFastaRecord(record.id, sampleName)
                record.description = record.id
                SeqIO.write(record, renamed, 'fasta')

        # subtyping
        subtype = getSubtype(folder)
        with open(out_results, 'a+') as results:
            results.write( str(Path(folder).name) + "\t" + "".join(subtype) + "\n" )
        # create subtype folder structure
        outSubtypeLoc = makeSubtypeFolderStruct(out_subtyped, subtype)
        os.makedirs(outSubtypeLoc, exist_ok=True)
        # copy fasta file
        shutil.copy(allFastaRenamed, outSubtypeLoc.joinpath(allFastaRenamed.name))
        
if __name__ == "__main__":
    run()