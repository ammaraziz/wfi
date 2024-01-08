# rule filter:
#     input:
#         single = INDIR / "{sample}.fastq.gz",
#     output:
#         filtered = OUTDIR / "qualtrim/{sample}.fastq",
#         status = OUTDIR / "status/filter_{sample}.txt"
#     params:
#         Fadapter = f"{ADAPTERS}/{ORG}_left.fa",
#         Radapter = f"{ADAPTERS}/{ORG}_right.fa"
#     threads: 2
#     message: "Trimming {wildcards.sample} reads."
#     log: OUTDIR / "logs/trim_{sample}.txt"
#     shell:"""
#         cutadapt {input.single}\
#         -j {threads} \
#         -g file:{params.Fadapter} \
#         -a file:{params.Radapter} \
#         -o {output.filtered} \
#         --report full &> {log}

#         touch {output.status}
#     """

rule trim:
    message: "Running fastp on {wildcards.sample}"
    input:
        single = INDIR / "{sample}.fastq.gz",
    output:
        filtered = OUTDIR / "qualtrim/{sample}.fastq.gz",
        html = OUTDIR / "qualtrimp" / "{sample}.html",
        json = OUTDIR / "qualtrimp" / "{sample}.json",
        status = OUTDIR / "status" / "filter_trim{sample}.txt",
    params:
        title = lambda wildcards: wildcards.sample,
    threads: 4
    log: OUTDIR / "logs/trimp_{sample}.txt"
    conda: "../envs/irma.yaml"
    shell:"""
    fastp \
    -i {input.r1} \
    -o {output.r1} \
    --json {output.json} \
    --html {output.html} \
    --report_title {params.title} \
    -w {threads} 2> {log}

    touch {output.status}
    """
