rule filter_std:
    input:
        faR1 = expand(config["input_dir"] + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
        faR2 = expand(config["input_dir"] + "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"])
    output:
        R1out = config["output_dir"] + "qualtrim/{sample}.R1.fastq",
        R2out = config["output_dir"] + "qualtrim/{sample}.R2.fastq",
        status = config["output_dir"] + "status/filter_{sample}.txt"
    params:
        Fadapter = f"{config['adapters']}/{org}_left.fa",
        Radapter = f"{config['adapters']}/{org}_right.fa"
    threads: 2
    message: "Filtering and trimming {wildcards.sample} reads."
    log: config["output_dir"] + "logs/trim_{sample}.txt"
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
        R1out = config["output_dir"] + "qualtrim/{sample}.R1.fastq",
        R2out = config["output_dir"] + "qualtrim/{sample}.R2.fastq"
    output:
        status = config["output_dir"] + "status/irma_{sample}.txt"
    params:
        folder = config["output_dir"] + "assemblies/{sample}/",
        run_module = lambda wildcards: irma_module,
    log: config["output_dir"] + "logs/irma_{sample}.txt"
    message: "IRMA is running for {wildcards.sample}"
    threads: 10
    shell:"""
        # run IRMA
        IRMA {params.run_module} {input.R1out} {input.R2out} {params.folder} 1>> {log}

        touch {output.status}
"""
