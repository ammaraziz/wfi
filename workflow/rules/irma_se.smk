rule irma:
    message: "IRMA is running for {wildcards.sample}"
    input:
        single = rules.time.output.filtered
    output:
        folder = directory(OUTDIR / "irma/{sample}/"),
        status = OUTDIR / "status/irma_{sample}.txt"
    params:
        run_module = lambda wildcards: IRMAMODULE,
    log: OUTDIR / "logs/irma_{sample}.txt"
    conda: "../envs/irma.yaml"
    threads: 10
    shell:"""
    IRMA {params.run_module} {input.single} {output.folder} 1> {log}
    touch {output.status}
    """