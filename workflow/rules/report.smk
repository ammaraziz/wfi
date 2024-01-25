rule renameSubtype:
    input:
        expand(config["output_dir"] + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        status = join(config["output_dir"] + "status/subtyping_complete.txt"),
    params:
        assemblies = config["output_dir"] + "irma/",
        output = config["output_dir"] + "assemblies/",
        org = ORG,
        subset = SEGTOKEEP
    shell:"""
    python scripts/subtyper.py \
    -i {params.assemblies:q} \
    -o {params.output:q} \
    -r {params.org} \
    --subset {params.subset}

    touch {output.status}
    """

rule SummaryReport:
    input:
        expand(config["output_dir"] + "status/irma_{sample}.txt", sample = SAMPLES)
    output:
        loc = join(config["output_dir"] + "status/plotting_complete.txt")
    params:
        ws = config["output_dir"] + "irma/",
        #loc_out = join(config["output_dir"] + "irma/"),
        loc_out = config["output_dir"] + "assemblies/",
        org = ORG
    shell:"""
    scripts/summaryReport.R \
    -i {params.ws:q} \
    -o {params.loc_out:q} \
    -r '{params.org}'

    touch {output.loc}
    """

rule stats:
    input:
        json = rules.mask_fasta.output.json,
    output: 
        json = WORKDIR / "irma" / "{sample}" / "json" / "{sample}.stats.json",
        status = WORKDIR / "status" / "stats_{sample}.txt",
    threads: 1
    conda: "../envs/irmakit.yaml"
    shell:"""
    python scripts/irmakit.py stats \
    --json {input.json} \
    --json-out {output.json}

    touch {output.status}
    """

rule plot:
    input:
        json = rules.stats.output.json,
        status = WORKDIR / "status" / "mask_{sample}.txt"
    output:
        html = WORKDIR / "irma" / "{sample}" / "{sample}.report.html",
        json = WORKDIR / "irma" / "{sample}" / "json" / "{sample}.plot.json",
        status = WORKDIR / "status" / "plot_{sample}.txt",
    params:
        sample_name = lambda wildcards: wildcards.sample,
        encoding = config['encoding']
    threads: 1
    conda: "../envs/irmakit.yaml"
    shell:"""
    python scripts/irmakit.py plot \
    --json {input.json} \
    --sample-name {params.sample_name} \
    --output-file {output.html} \
    --encoding {params.encoding} \
    --json-out {output.json}
    
    touch {output.status}
    """
