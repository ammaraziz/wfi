rule irma:
    message: "IRMA is running for {wildcards.sample}"
    input:
        single = rules.trim.output.filtered
    output:
        status = WORKDIR / "status" / "irma_{sample}.txt"
    params:
        irma_dir = WORKDIR / "irma" / "{sample}/",
        run_module = lambda wildcards: IRMAMODULE,
    log: WORKDIR / "logs/irma_{sample}.txt"
    conda: "../envs/irma.yaml"
    threads: 10
    shell:"""
    mkdir -p {params.irma_dir}

    IRMA {params.run_module} {input.single} {params.irma_dir} 1> {log}
    
    touch {output.status}
    """