rule irma:
    message: "IRMA is running for {wildcards.sample}"
    input:
        r1 = rules.trim.output.r1,
        r2 = rules.trim.output.r2,
    output:
        status = WORKDIR / "status/irma_{sample}.txt",
    params:
        irma_dir = lambda wildcards: WORKDIR / "irma" / wildcards.sample,
        run_module = IRMAMODULE,
    conda: "../envs/irma.yaml"
    log: WORKDIR / "logs/irma_{sample}.txt"
    threads: 10
    shell:"""
    
    IRMA {params.run_module} {input.r1} {input.r2} {params.irma_dir} 1> {log}

    touch {output.status}
    """
