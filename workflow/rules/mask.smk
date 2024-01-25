rule mask_fasta:
    input:
        irmastatus = rules.irma.output.status,
        json = rules.irma_parse.output.json
    output:
        status = WORKDIR / "status" / "mask_{sample}.txt",
        json = WORKDIR / "irma" / "{sample}" / "json" / "{sample}.mask.json",
    params:
        sample_name = lambda w: w.sample,
        outdir = lambda wildcards: WORKDIR / "irma" / wildcards.sample / "masked/",
        min_cov = config['maskdepth'] 
    threads: 1
    conda: "../envs/irmakit.yaml"
    params:
        IRMA = WORKDIR
    shell:"""
    python scripts/irmakit.py maskfasta \
    --json {input.json} \
    --sample-name {params.sample_name} \
    --output-dir {params.outdir} \
    --min-cov {params.min_cov} \
    --json-out {output.json}

    touch {output.status}
    """