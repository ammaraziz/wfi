#############################################################################
#                                                                           #
#                           wfi  (WHO-FLU-IRMA)                             #
#         pipeline for the assembly of illumina short read data             #
#                                                                           #
#                          Created by Ammar Aziz                            #
#                                                                           #
#############################################################################


#############################################################################
#                           DO NOT TOUCH ANYTHING                           #
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
#trimmomatic = config["trimmomatic"]
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
            seg_to_keep = "{HA,NA,MP}"
            #seg_to_keep = 'subset'
        elif subset is False:
            seg_to_keep = "{HA,NA,MP,NS,NP,PA,PB1,PB2}"
            #seg_to_keep = 'all'
        else:
            raise ValueError("Check config file for 'subset' param. If unsure set to: False")

    elif org == 'RSV':
        #print('Organism {}'.format(org), file = sys.stdout)
        mode = 'RSV'
        seg_to_keep = "rsv_"
        #seg_to_keep = "rsv"

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
            segment = sub('^[A|B]_', '', record.description)
            # for g in listSeg:
            #     if g in record.description:
            #         segment = g
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
        R2out = workspace + "qualtrim/{sample}.R2.paired.fastq"
        #R1out_unpaired = workspace + "qualtrim/{sample}.R1.unpaired.fastq",
        #R2out_unpaired = workspace + "qualtrim/{sample}.R2.unpaired.fastq"
    params:
       #trimmo = trimmomatic,
       #cutadapt = config["cutadapt"]
       Fadapter = f"bin/adapters/{org}_f.fa",
       Radapter = f"bin/adapters/{org}_r.fa"
    threads: 2
    message: "Filtering and trimming {input.faR1} reads."
    log: workspace + "logs/trim_{sample}.txt"
    shell:"""
      #cutadapt
      cutadapt {input.faR1} {input.faR2} -j {threads} -g file:{params.Fadapter} -A file:{params.Radapter} -o {output.R1out} -p {output.R2out} --report full 1> {log}
    """

#Assembly using IRMA PE mode.
checkpoint irma:
    input:
        R1out = workspace + "qualtrim/{sample}.R1.paired.fastq",
        R2out = workspace + "qualtrim/{sample}.R2.paired.fastq"
    output:
        contigs = workspace + "assemblies/{sample}/contigs.fasta",
        status = temp(workspace + "assemblies/{sample}/irma_status.txt")
    params:
        sample_name = "{sample}",
        afolder = workspace + "assemblies/",
        folder = workspace + "assemblies/{sample}/",
        segs = lambda widlcards: seg_to_keep,
        run_mode = lambda wildcards: mode,
        vcf_loc = workspace + 'vcf/' + "{sample}/"
    log: workspace + "logs/irma_{sample}.txt"
    message: "IRMA is running for {input.R1out}"
    threads: 10
    shell:"""
        # run IRMA
        IRMA {params.run_mode} {input.R1out} {input.R2out} {params.sample_name} 1>> {log}
        # move output to folder
        mv $PWD/{params.sample_name} {params.afolder}

        myarray=(`find ./ -maxdepth 1 -name "*.fasta"`)
        if [ ${{#myarray[@]}} -gt 0 ]
        echo ${{myarray[@]}}
        then
            cat {params.folder}*.fasta 1> {output.contigs} 2>> {log}
            #python tools/geneMover.py {params.folder} {output.contigs} {params.segs} 2>> {log}
            cat {params.folder}amended_consensus/*.fa 1> {params.folder}amended_consensus/amended.contigs.fasta 2>>{log}
            mkdir -p {params.vcf_loc} && cp {params.folder}*.vcf {params.vcf_loc}
            echo irma produced stuff this is temp ignore > {output.status}
        else
            echo temp file ignore > {output.contigs}
            touch temp file ignore > {output.status}
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
        loc = workspace + "assemblies/rename/{sample}_"
    params:
        ws = workspace + "assemblies/",
        org = org
    shell:"""
        ./tools/summaryReport.R -i {params.ws:q} -o {output.loc:q} -r '{params.org}'
    """

# TODO

# Fix issue with FLU, wifi not moving fasta and renaming.