rule trim:
    message: "Filtering and trimming {wildcards.sample} reads."
    input:
        # r1 = INDIR / "{sample}_L001_R1_001.fastq.gz",
        # r2 = INDIR / "{sample}_L001_R2_001.fastq.gz",
        r1 = get_input_r1,
        r2 = get_input_r2,
    output:
        r1 = OUTDIR / "qualtrimp" / "{sample}.R1.fastq.gz",
        r2 = OUTDIR / "qualtrimp" / "{sample}.R2.fastq.gz",
        html = OUTDIR / "qualtrimp" / "{sample}.html",
        json = OUTDIR / "qualtrimp" / "{sample}.json",
        status = OUTDIR / "status" / "fastp_{sample}.txt",
    params:
        title = lambda wildcards: wildcards.sample,
    threads: 4
    log: OUTDIR / "logs/trimp_{sample}.txt"
    conda: "../envs/trim.yaml"
    shell:"""
    fastp \
    -i {input.r1} \
    -I {input.r2} \
    -o {output.r1} \
    -O {output.r2} \
    --detect_adapter_for_pe \
    --json {output.json} \
    --html {output.html} \
    --report_title {params.title} \
    -w {threads} 2> {log}

    touch {output.status}
    """
