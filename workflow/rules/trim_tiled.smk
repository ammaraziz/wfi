rule trim:
    message: "Trimming tiled amplicons for {wildcards.sample}"
    input:
        r1 = INDIR / "{sample}_L001_R1_001.fastq.gz",
        r2 = INDIR / "{sample}_L001_R2_001.fastq.gz",
    output:
        r1 = OUTDIR / "trim_tile" / "{sample}.R1.fastq.gz",
        r2 = OUTDIR / "trim_tile" / "{sample}.R2.fastq.gz",
        status = OUTDIR / "status" / "trim_{sample}.txt",
        summary = OUTDIR / "trim_tile" / "{sample}.ptrimmer.summary.txt",
    params:
        adapters = f"resources/adapters/{ORG}_tiled.fa",
        minqual = 15,
        kmer = 8,
        mismatch = 1,
        r1tmp = OUTDIR / "trim_tile" / "{sample}.R1.fastq",
        r2tmp = OUTDIR / "trim_tile" / "{sample}.R2.fastq",
    threads: 1
    conda: "../envs/tiled.yaml"
    log: OUTDIR / "logs" / "tiled_{sample}.txt"
    shell:"""
    ptrimmer \
    --keep \
    --seqtype pair \
    --ampfile {params.adapters} \
    --summary {output.summary} \
    --minqual {params.minqual} \
    --kmer {params.kmer} \
    --mismatch {params.mismatch} \
    --read1 {input.r1} \
    --trim1 {params.r1tmp} \
    --read2 {input.r2} \
    --trim2 {params.r2tmp} &> {log}

    gzip -c {params.r1tmp} > {output.r1}
    gzip -c {params.r2tmp} > {output.r2}

    rm {params.r1tmp}
    rm {params.r2tmp}
    
    touch {output.status}
    """