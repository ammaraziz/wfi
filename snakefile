#############################################################################
#                                                                           #
#                           wfi  (WHO-FLU-IRMA)                             #
#   pipeline for the assembly of illumina short read data                   #
#                                                                           #
#   Created by Miguel Lopez with heavy modification by Ammar Aziz,          #
#                   helper script by Uma. #                                 #
#                                                                           #
#############################################################################


#############################################################################
#                    DO NOT TOUCH ANYTHING BELOW THIS LINE                  #
#############################################################################


import subprocess, sys, os, glob, shutil
from re import sub
from os.path import join #needed?
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# export IRMA into $PATH of linux
irma_path = "bin/flu-amd/"
os.environ["PATH"] += os.pathsep + os.pathsep.join([irma_path])


# environment variables
configfile: "wfi_config.yaml"
IFQ = config["input_dir"]
workspace = config["output_dir"]
trimmomatic = config["trimmomatic"]
org = config["organism"].upper()
subset = config["subset"] #keep only ha, na, mp
secondary_assembly = config["second_assembly"]

# set organism and gene segements (influenza) to keep
if org in ['FLU', 'RSV']:
    if org == 'FLU':

        # primary/secondary assembly settings
        if secondary_assembly is True:
            mode = 'FLU-secondary'
        elif secondary_assembly is False:
            mode = 'FLU'
        else:
            raise ValueError("Assembly mode unknown. Check config file for 'second_assembly', options MUST be either True or False. (exactly)")

        # gene segment settings
        if subset is True:
            #seg_to_keep = "{HA,NA,MP}"
            seg_to_keep = 'subset'
        elif subset is False:
            #seg_to_keep = "{HA,NA,MP,NS,NP,PA,PB1,PB2}"
            seg_to_keep = 'all'
        else:
            raise ValueError("Check config file for 'subset' param. If unsure set to: False")

    elif org == 'RSV':
        #print('Organism {}'.format(org), file = sys.stdout)
        mode = 'RSV'
        #seg_to_keep = "rsv_"
        seg_to_keep = "rsv"

    else:
        raise ValueError("Check config file for 'organism' setting. Options are: FLU or RSV")

## Functions -------------------------------------------------------------------

def fixNames(fafile):

    if org == 'FLU':
        res = ""
        segment = ""
        listSeg = ["HA","MP","NA","NP","NS","PA","PB1","PB2"]
        sampleName = fafile.split("/")[-2].split("_")[0]
        for index, record in enumerate(SeqIO.parse(fafile, "fasta")):
            for g in listSeg:
                if g in record.description:
                    segment = g
            res += ">" + sampleName + "_" + segment + "\n" + str(record.seq) + "\n"
        return(res)

    elif org == 'RSV':
        res = ""
        sampleName = fafile.split("/")[-2].split("_")[0]
        for index, record in enumerate(SeqIO.parse(fafile, "fasta")):
            res = ">" + sampleName + "\n" + str(record.seq) + "\n"
        return(res)

    else:
        raise ValueException("Org not found. Error: fixNames")


## Rules ------------------------------------------------------------------------

#PAIRED READS
SAMPLES, PAIR = glob_wildcards(IFQ + "/{sample}_L001_{pair}_001.fastq.gz")

rule all:
    input:
        expand(workspace + "qualtrim/{sample}.R1.paired.fastq", sample = SAMPLES),
        expand(workspace + "qualtrim/{sample}.R2.paired.fastq", sample = SAMPLES),
        expand(workspace + "assemblies/{sample}/contigs.fasta", sample = SAMPLES),
        expand(workspace + "assemblies/rename/{sample}.fasta", sample = SAMPLES),
        expand(workspace + "assemblies/rename/{sample}.txt", sample = SAMPLES),
        expand(workspace + "assemblies/{sample}/irma_status.txt", sample = SAMPLES),
        expand(workspace + "logs/run_report.pdf")

#QUALITY FILTER
rule filter:
    input:
        faR1 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
        faR2 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"])
    output:
        R1out = workspace + "qualtrim/{sample}.R1.paired.fastq",
        R2out = workspace + "qualtrim/{sample}.R2.paired.fastq",
        R1out_unpaired = workspace + "qualtrim/{sample}.R1.unpaired.fastq",
        R2out_unpaired = workspace + "qualtrim/{sample}.R2.unpaired.fastq"
    params:
       trimmo = trimmomatic
    threads: 1
    message: "Filtering and trimming {input.faR1} reads."
    log: workspace + "logs/trimmomatic_{sample}.txt"
    shell:"""
      java -jar {params.trimmo}/trimmomatic-0.39.jar PE -threads {threads} -phred33 {input.faR1} {input.faR2} {output.R1out} {output.R1out_unpaired} {output.R2out} {output.R2out_unpaired} ILLUMINACLIP:{params.trimmo}/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:75 HEADCROP:10 TRAILING:5 2> {log}
    """

#Assembly using IRMA PE mode.
checkpoint irma:
    input:
        R1out = workspace + "qualtrim/{sample}.R1.paired.fastq",
        R2out = workspace + "qualtrim/{sample}.R2.paired.fastq"
    output:
        contigs = workspace + "assemblies/{sample}/contigs.fasta",
        status = workspace + "assemblies/{sample}/irma_status.txt"
    params:
        sample_name = "{sample}",
        afolder = workspace + "assemblies/",
        folder = workspace + "assemblies/{sample}/",
        segs = lambda widlcards: seg_to_keep,
        run_mode = lambda wildcards: mode
    log: workspace + "logs/irma_{sample}.txt"
    message: "IRMA is running for {input.R1out}"
    threads: 10
    shell:"""
        IRMA {params.run_mode} {input.R1out} {input.R2out} {params.sample_name} 1>> {log}
        mv $PWD/{params.sample_name} {params.afolder}
        if [ -s {params.folder}*.fasta ]
        then
            cat {params.folder}*{params.segs}*.fasta 1> {output.contigs} 2>> {log}
            # python tools/geneMover.py {params.folder} {output.contigs} {params.segs} 2>> {log}
            cat {params.folder}amended_consensus/*.fa 1> {params.folder}amended_consensus/amended.contigs.fasta 2>>{log}
            touch {output.status}
        else
            touch {output.contigs}
            touch {output.status}
        fi
    """

rule rename_fasta:
    input:
        fasta = workspace + "assemblies/{sample}/contigs.fasta"
    output:
        fasta = workspace + "assemblies/rename/{sample}.fasta"
    run:
        newFasta  =  fixNames(input.fasta)
        with open(output.fasta, "w") as fo:
            fo.write(newFasta)


# sort into subtypes
rule subTypeSort:
    input:
        fasta = workspace + "assemblies/rename/{sample}.fasta"
    output:
        tmp = workspace + "assemblies/rename/{sample}.txt"
    params:
        ws = workspace,
        fasta = workspace + "assemblies/{sample}/contigs.fasta",
        org = org
    run:
        if params.org == 'FLU':
            shell: """
                python ./tools/groupSubtypeIRMA.py {params.fasta} {params.ws}
                echo temp file ignore > {output.tmp}
            """

        elif params.org == 'RSV':
            #tmp function; replace with python code to sort rsv subtypes into individual
            with open(output.tmp, "w") as file_open:
                file_open.write("tmp file ignore")

        else:
            shell: """
            touch {output.tmp}
            """

rule SummaryReport:
    input:
        expand(workspace + "assemblies/{sample}/irma_status.txt", sample = SAMPLES)
    output:
        pdf = workspace + "logs/run_report.pdf"
    params:
        ws = workspace,
        org = org
        #table = expand(workspace + "assemblies/{sample}/tables/READ_COUNTS.txt", sample = SAMPLES)
    shell:"""
        ./tools/irma_summary.R -i {params.ws:q} -o {output.pdf:q} -r '{params.org}'
    """

# TODO

#wfi to do list:
#   Fix issue with H3 275Y position (it should be 274)
    # similar issue with Bvic (do not call 275Y )
#•  Fix output of subset so B virus only includes HA/NA (not MP)
#   Fix issue with some genes missing causing an error and snakemake pipeline to fail
#•  Create a summary file to act as a report for quick diagnostics, includes:
#   o   Average depth for each gene
#   o   Presence of ambiguous bases in amended files
#   o   Suspected mixture if observed read stats are suspicious
#