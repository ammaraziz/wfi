rule trim:
    message: "Running fastp on {wildcards.sample}"
    input:
        single = INDIR / "{sample}.fastq.gz",
    output:
        filtered = OUTDIR / "qualtrim" / "{sample}.fastq.gz",
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
