rule trim:
    message: "Running fastp on {wildcards.sample}"
    input:
        single = INDIR / "{sample}.fastq.gz",
    output:
        filtered = WORKDIR / "qc" / "{sample}.fastq.gz",
        html = WORKDIR / "qc" / "{sample}.html",
        json = WORKDIR / "qc" / "{sample}.json",
        status = WORKDIR / "status" / "filter_trim{sample}.txt",
    params:
        title = lambda wildcards: wildcards.sample,
    threads: 4
    log: WORKDIR / "logs" / "trimp_{sample}.txt"
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
