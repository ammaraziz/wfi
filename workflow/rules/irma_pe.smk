rule irma:
    message: "IRMA is running for {wildcards.sample}"
    input:
        r1 = rules.trim.output.r1,
        r2 = rules.trim.output.r2,
    output:
        status = OUTDIR / "status/irma_{sample}.txt",
        folder = directory(OUTDIR / "irma/{sample}/"),
    params:
        run_module = IRMAMODULE,
    conda: "../envs/irma.yaml"
    log: OUTDIR / "logs/irma_{sample}.txt"
    threads: 10
    shell:"""
    IRMA {params.run_module} {input.r1} {input.r2} {output.folder} 1> {log}

    touch {output.status}
    """
