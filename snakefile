#############################################################################
#                                                                           #
#                           wfi  (WHO-FLU-IRMA)                             #
#                      One Big IRMA Wrapper in Snakemake                    #
#                                                                           #
#                          Created by Ammar Aziz                            #
#                                                                           #
#############################################################################


#############################################################################
#                           DO NOT TOUCH ANYTHING                           #
#############################################################################

version = "0.3.5"

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
trim_prog = config["trim_prog"]
trim_org = config["trim_org"]

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

 ## Message ------------------------------------------------------------------------
onstart:
    print("Run mode: " + run_mode)
    print("Sequence Technology: " + seq_technology)
    print("Organism: " + org)
    print("IRMA Module: " + irma_module)
    print("Secondary Assembly: " + str(secondary_assembly))
    print("Trimming using: " + str(trim_prog))
    print("\n")

onsuccess:
    print("wfi has successfully completed!")

onerror:
    print("oops wfi has run into an issue. Look above, If rule SummaryReport failed there is no need to worry!")

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
if run_mode not in ['paired', 'single']:
    sys.exit("Configuration incorrect, check 'run_mode' it must be: paired or single")
if run_mode == 'paired':
    rule_mode = [
    expand(workspace + "qualtrim/{sample}.R1.fastq", sample = SAMPLES),
    expand(workspace + "qualtrim/{sample}.R2.fastq", sample = SAMPLES)
    ]
if run_mode == 'single':
    rule_mode = [
    expand(workspace + "qualtrim/{sample}.fastq", sample = SAMPLES)
    ]

## Trimming ---------------------------------------------------------------------
if trim_prog not in ['standard', 'tile']:
    sys.exit("Configuration incorrect, check 'trim_prog' it must be: standard or tile")
if trim_org not in ['h1', 'h3']:
    sys.exit("Configuration incorrect, check 'trim_org' it must be: h1 or h3. bvic is not supported")


## Rules ------------------------------------------------------------------------
rule all:
    input:
        # singe or paired
        rule_mode,
        # status
        expand(workspace + "status/filter_{sample}.txt", sample = SAMPLES),
        expand(workspace + "status/irma_{sample}.txt", sample = SAMPLES),
        join(workspace + "status/subtyping_complete.txt"),
        join(workspace + "status/plotting_complete.txt")



if run_mode == 'paired':
    # Filter standard
    if trim_prog == 'standard':
        rule filter_std:
            input:
                faR1 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
                faR2 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"])
            output:
                R1out = workspace + "qualtrim/{sample}.R1.fastq",
                R2out = workspace + "qualtrim/{sample}.R2.fastq",
                status = workspace + "status/filter_{sample}.txt"
            params:
               Fadapter = f"bin/adapters/{org}_left.fa",
               Radapter = f"bin/adapters/{org}_right.fa"
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
    # Filter tile
    if trim_prog == 'tile':
        rule filter_tile:
            input:
                faR1 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
                faR2 = expand(IFQ + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"])
            output:
                R1_out_left = workspace + "qualtrim/{sample}_left.R1.fastq",
                R2_out_left = workspace + "qualtrim/{sample}_left.R2.fastq",
                R1_out_right = workspace + "qualtrim/{sample}.R1.fastq",
                R2_out_right = workspace + "qualtrim/{sample}.R2.fastq",
                status = workspace + "status/filter_{sample}.txt"
            params:
               Fadapter = f"bin/adapters/{org}_tile_left_{trim_org}.fa",
               Radapter = f"bin/adapters/{org}_tile_right_{trim_org}.fa",
               stats1 = workspace + "logs/trimStats1_{sample}.txt",
               stats2 = workspace + "logs/trimStats2_{sample}.txt",
               refstats1 = workspace + "logs/trimRefStats1_{sample}.txt",
               refstats2 = workspace + "logs/trimRefStats2_{sample}.txt",
               k = 9,
               mink = 3,
               restrict = 30,
               hdist = 1
            threads: 2
            message: "Filtering and trimming {input.faR1} reads."
            log: workspace + "logs/trim_{sample}.txt"
            shell:"""
              bbduk.sh in={input.faR1} \
              in2={input.faR2} \
              ktrim=l \
              mm=f \
              hdist={params.hdist} \
              rcomp=t \
              ref={params.Fadapter} \
              ordered=t \
              minlen=0 \
              minlength=0 \
              trimq=0 \
              k={params.k} \
              mink={params.mink} \
              threads={threads} \
              restrictleft={params.restrict} \
              out={output.R1_out_left} \
              out2={output.R2_out_left} \
              stats={params.stats1} \
              statscolumns=5 2>> {log}
              
              bbduk.sh in={output.R1_out_left} \
              in2={output.R2_out_left} \
              ktrim=r \
              mm=f \
              hdist={params.hdist} \
              rcomp=t \
              ref={params.Radapter} \
              ordered=t \
              minlen=0 \
              minlength=0 \
              trimq=0 \
              k={params.k} \
              mink={params.mink} \
              threads={threads} \
              restrictright={params.restrict} \
              out={output.R1_out_right} \
              out2={output.R2_out_right} \
              stats={params.stats2} \
              statscolumns=5 2>> {log}

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
            folder = workspace + "assemblies/{sample}/",
            run_module = lambda wildcards: irma_module,
        log: workspace + "logs/irma_{sample}.txt"
        message: "IRMA is running for {input.R1out}"
        threads: 10
        shell:"""
            # run IRMA
            IRMA {params.run_module} {input.R1out} {input.R2out} {params.folder} 1>> {log}

            touch {output.status}
    """

elif run_mode == 'single':
    rule filter_std:
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
            folder = workspace + "assemblies/{sample}/",
            run_module = lambda wildcards: irma_module,
        log: workspace + "logs/irma_{sample}.txt"
        message: "IRMA is running for {input.single}"
        threads: 10
        shell:"""
            IRMA {params.run_module} {input.single} {params.folder} 1>> {log}
            touch > {output.status}
        """
else: 
    sys.exit("Something went wrong with the filter+irma command. Check output")

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
