rule irma_parse:
    input:
        irmastatus = rules.irma.output.status,
    output:
        json = WORKDIR / "irma" / "{sample}" / "json" / "{sample}.json",
        status = WORKDIR / "status" / "parse_{sample}.txt"
    params:
        sample_name = lambda widlcards: widlcards.sample,
        irma_dir = rules.irma.params.irma_dir
    conda: "../envs/irmakit.yaml"
    threads: 1
    shell:"""
    python scripts/irmakit.py parse \
    --irma-dir {params.irma_dir} \
    --json {output.json} \
    --sample-name {params.sample_name}

    touch {output.status}
    """

rule irma_stats:
    input:
        json = rules.irma_parse.output.json
    output:
        status = WORKDIR / "status" / "stats_{sample}.txt"
    conda: "..envs/irmakit.yaml"
    threads: 1
    shell:"""
    python scripts/irmakit.py stats \
    --json {input.json} \

    touch {output.status}
    """