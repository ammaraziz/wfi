rule filter:
    input:
        single = expand(config["input_dir"] + "{{sample}}.fastq.gz"),
    output:
        filtered = config["output_dir"] + "qualtrim/{sample}.fastq",
        status = config["output_dir"] + "status/filter_{sample}.txt"
    params:
        Fadapter = f"{config['adapters']}/{ORG}_left.fa",
        Radapter = f"{config['adapters']}/{ORG}_right.fa"
    threads: 2
    message: "Trimming {wildcards.sample} reads."
    log: config["output_dir"] + "logs/trim_{sample}.txt"
    shell:"""
        cutadapt {input.single}\
        -j {threads} \
        -g file:{params.Fadapter} \
        -a file:{params.Radapter} \
        -o {output.filtered} \
        --report full &> {log}

        touch {output.status}
    """

rule irma_single:
    input:
        single = config["output_dir"] + "qualtrim/{sample}.fastq"
    output:
        status = config["output_dir"] + "status/irma_{sample}.txt"
    params:
        folder = config["output_dir"] + "irma/{sample}/",
        run_module = lambda wildcards: IRMAMODULE,
    log: config["output_dir"] + "logs/irma_{sample}.txt"
    message: "IRMA is running for {wildcards.sample}"
    threads: 10
    shell:"""
        IRMA {params.run_module} {input.single} {params.folder} 1> {log}
        touch {output.status}
    """