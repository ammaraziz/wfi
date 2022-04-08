#############################################################################
#                                                                           #
#                           wfi  (WHO-FLU-IRMA)                             #
#                   pipeline for the assembly of illumina/ont               #
#                                                                           #
#                          Created by Ammar Aziz                            #
#                                                                           #
#############################################################################


#############################################################################
#                           DO NOT TOUCH ANYTHING                           #
#############################################################################

version = "0.3.3"

import subprocess, sys, os, glob, shutil
from time import sleep
from re import sub
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import snakemake

# export IRMA into $PATH of linux
irma_path = "bin/flu-amd/"
os.environ["PATH"] += os.pathsep + os.pathsep.join([irma_path])


# global variables
configfile: "wfi_config.yaml"
IFQ = config["input_dir"]
workspace = config["output_dir"]
org = config["organism"].upper()
seq_technology = config["technology"].lower()
secondary_assembly = config["secondary_assembly"]
subset = config["subset"]

# set organism and gene segements (influenza) to keep

if org == 'FLU':
    if seq_technology == 'illumina':
        if secondary_assembly is True:
            irma_module = 'FLU-secondary'
        else:
            irma_module = 'FLU'
    elif seq_technology == 'ont' and secondary_assembly == False:
        irma_module = 'FLU-minion'
    elif seq_technology == 'pgm' and secondary_assembly == False:
        irma_module = 'FLU-pgm'
    else:
        sys.exit(
            f'Assembly module unknown. Check config for options:\n'
            f'organism {org}\n' 
            f'secondary_assembly {secondary_assembly}\n'
            f'techology {seq_technology}'
            )

    # gene segment settings
    if subset is True:
        seg_to_keep = "subset"
    elif subset is False:
        seg_to_keep = "all"
    else:
        sys.exit("Check config file for 'subset' param. If unsure set to: False")

elif org == 'RSV':
    seg_to_keep = "rsv"
    if seq_technology == 'illumina': 
        if secondary_assembly is True:
            irma_module = 'RSV-secondary'
        else:
            irma_module = 'RSV'
    elif seq_technology and secondary_assembly is False:
        irma_module = 'RSV-minion'
    else:
        sys.exit(
            f'Assembly module unknown. Check config for options:\n'
            f'organism {org}\n' 
            f'secondary_assembly {secondary_assembly}\n'
            f'techology {seq_technology}'
            )

else:
    raise ValueError("Check config file for 'organism' setting. Options are: FLU or RSV")

onsuccess:
    print("wfi has successfully completed!")

onerror:
    print("oops wfi has run into an issue. Look above, If rule SummaryReport failed there is no need to worry!")

## Functions -------------------------------------------------------------------

def fixNames(fafile, name, org):
    sample_name, sample_number = name.split("_")

    if org == 'FLU':
        res = ""
        listSeg = {"HA":"4","MP":"7","NA":"6","NP":"5","NS":"8","PA":"3","PB1":"2","PB2":"1"}
        for index, record in enumerate(SeqIO.parse(fafile, "fasta")):
            gene = record.id.split("_")[1]
            res += ">" + sample_name + "." + listSeg[gene] + "\n" + str(record.seq) + "\n"
        return(res)

    elif org == 'RSV':
        res = ""
        for index, record in enumerate(SeqIO.parse(fafile, "fasta")):
            res = ">" + sample_name + "\n" + str(record.seq) + "\n"
        return(res)

    else:
        sys.exit("Org not found. Error: fixNames")

## Sequencing technology ------------------------------------------------------------------------
if seq_technology == 'illumina':
    SAMPLE_NAME, SAMPLE_NUMBER, lane_number, PAIR = glob_wildcards(IFQ + "/{sample_name}_{sample_number}_L{lane_number}_{pair}_001.fastq.gz")
    SAMPLES = [i + "_" + x for i, x in zip(SAMPLE_NAME, SAMPLE_NUMBER)]
    run_mode = 'paired'
elif seq_technology == 'ont':
    SAMPLE_NAME = glob_wildcards(IFQ + config['pattern_ont'])[0]
    SAMPLES = SAMPLE_NAME
    run_mode = 'single'
elif seq_technology == 'pgm':
    SAMPLE_NAME = glob_wildcards(IFQ + config['pattern_pgm'])[0]
    SAMPLES = SAMPLE_NAME
    run_mode = 'single'

## Mode ------------------------------------------------------------------------

if run_mode == 'paired':
    rule_mode = [expand(workspace + "qualtrim/{sample}.R1.fastq", sample = SAMPLES),
    expand(workspace + "qualtrim/{sample}.R2.fastq", sample = SAMPLES)]
elif run_mode == 'single':
    rule_mode = [expand(workspace + "qualtrim/{sample}.fastq", sample = SAMPLES)]

## Message ------------------------------------------------------------------------
print("Run mode: " + run_mode)
print("Sequence Technology: " + seq_technology)
print("Organism: " + org)
print("IRMA Module: " + irma_module)
print("Secondary Assembly: " + str(secondary_assembly))
print("\n")
## Rules ------------------------------------------------------------------------

rule all:
    input:
        # filter
        rule_mode,
        # status
        expand(workspace + "status/filter_{sample}.txt", sample = SAMPLES),
        expand(workspace + "status/irma_{sample}.txt", sample = SAMPLES),
        join(workspace + "status/subtyping_complete.txt"),
        join(workspace + "status/plotting_complete.txt")


if run_mode == 'paired':
    # Filter
    rule filter_paired:
        input:
            faR1 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
            faR2 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"])
        output:
            R1out = workspace + "qualtrim/{sample}.R1.fastq",
            R2out = workspace + "qualtrim/{sample}.R2.fastq",
            status = workspace + "status/filter_{sample}.txt"
        params:
           Fadapter = f"bin/adapters/{org}_f.fa",
           Radapter = f"bin/adapters/{org}_r.fa"
        threads: 2
        message: "Filtering and trimming {input.faR1} reads."
        log: workspace + "logs/trim_{sample}.txt"
        shell:"""
          cutadapt {input.faR1} {input.faR2} \
          -j {threads} \
          -g file:{params.Fadapter} \
          -A file:{params.Radapter} \
          -o {output.R1out} \
          -p {output.R2out} \
          --report full 1> {log}

          touch {output.status}
        """

    # Assembly
    rule irma_paired:
        input:
            R1out = workspace + "qualtrim/{sample}.R1.fastq",
            R2out = workspace + "qualtrim/{sample}.R2.fastq"
        output:
            status = workspace + "status/irma_{sample}.txt"
        params:
            sample_name = "{sample}",
            afolder = workspace + "assemblies/",
            folder = workspace + "assemblies/{sample}/",
            run_module = lambda wildcards: irma_module,
            vcf_loc = workspace + 'vcf/' + "{sample}/"
        log: workspace + "logs/irma_{sample}.txt"
        message: "IRMA is running for {input.R1out}"
        threads: 10
        shell:"""
            # run IRMA
            IRMA {params.run_module} {input.R1out} {input.R2out} {params.sample_name} 1>> {log}
            
            # move output to folder
            mv $PWD/{params.sample_name} {params.afolder}

            touch {output.status}
    """

elif run_mode == 'single':
    rule filter_single:
        input:
            single = expand(IFQ + "{{sample}}.fastq.gz"),
        output:
            filtered = workspace + "qualtrim/{sample}.fastq",
            status = workspace + "status/filter_{sample}.txt"
        params:
           Fadapter = f"bin/adapters/{org}_f.fa",
           Radapter = f"bin/adapters/{org}_r.fa"
        threads: 2
        message: "Trimming {input.single} reads."
        log: workspace + "logs/trim_{sample}.txt"
        shell:"""
          cutadapt {input.single}\
          -j {threads} \
          -g file:{params.Fadapter} \
          -a file:{params.Radapter} \
          -o {output.filtered} \
          --report full 1> {log}

          touch {output.status}
        """

    rule irma_single:
        input:
            single = workspace + "qualtrim/{sample}.fastq"
        output:
            status = workspace + "status/irma_{sample}.txt"
        params:
            sample_name = "{sample}",
            afolder = workspace + "assemblies/",
            folder = workspace + "assemblies/{sample}/",
            run_module = lambda wildcards: irma_module,
            vcf_loc = workspace + 'vcf/' + "{sample}/"
        log: workspace + "logs/irma_{sample}.txt"
        message: "IRMA is running for {input.single}"
        threads: 10
        shell:"""
            IRMA {params.run_module} {input.single} {params.sample_name} 1>> {log}
            mv $PWD/{params.sample_name} {params.afolder}
            touch > {output.status}
        """
else: 
    sys.exit("Something went wrong with the filter+irma command")

# rule rename_fasta:
#     input:
#         fasta = workspace + "assemblies/{sample}/contigs.fasta",
#     output:
#         fasta = workspace + "assemblies/rename/{sample}.fasta"
#     params:
#         org = org,
#         sample_name = "{sample}"
#     run:
#         newFasta = fixNames(fafile = input.fasta, name = params.sample_name, org = params.org)
#         with open(output.fasta, "w") as fo:
#             fo.write(newFasta)


# rename and sort into subtypes
rule renameSubtype:
    input:
        expand(workspace + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        status = join(workspace + "status/subtyping_complete.txt")
    params:
        assemblies = workspace + "assemblies/",
        org = org,
        subset = seg_to_keep
    shell:"""
    python tools/subtyper.py \
    -i {params.assemblies} \
    -o {params.assemblies} \
    -r {params.org} \
    --subset {params.subset:q}

    touch {output.status}
    """

rule SummaryReport:
    input:
        expand(workspace + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        loc = join(workspace + "status/plotting_complete.txt")
    params:
        ws = workspace + "assemblies/",
        loc_out = join(workspace + "assemblies/"),
        org = org
    shell:"""
        ./tools/summaryReport.R \
        -i {params.ws:q} \
        -o {params.loc_out:q} \
        -r '{params.org}'
        
        touch {output.loc}
    """
