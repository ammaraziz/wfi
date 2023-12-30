rule renameSubtype:
    input:
        expand(config["output_dir"] + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        status = join(config["output_dir"] + "status/subtyping_complete.txt")
    params:
        assemblies = config["output_dir"] + "assemblies/",
        org = ORG,
        subset = SEGTOKEEP
    shell:"""
    python scripts/subtyper.py \
    -i {params.assemblies} \
    -o {params.assemblies} \
    -r {params.org} \
    --subset {params.subset:q}

    touch {output.status}
    """

rule SummaryReport:
    input:
        expand(config["output_dir"] + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        loc = join(config["output_dir"] + "status/plotting_complete.txt")
    params:
        ws = config["output_dir"] + "assemblies/",
        loc_out = join(config["output_dir"] + "assemblies/"),
        org = ORG
    shell:"""
        scripts/summaryReport.R \
        -i {params.ws:q} \
        -o {params.loc_out:q} \
        -r '{params.org}'

        touch {output.loc}
    """