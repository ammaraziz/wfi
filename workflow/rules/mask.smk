rule createBedForMask:
    input:
        irmaDir = directory(OUTDIR / "irma" / "{sample}/")
    output:
        maskbed = OUTDIR / "mask" / "{sample}" / "{sample}.mask.bed"
    run:
        irma_files = getIRMAFiles(input)
        

rule maskfasta:
    input:
        ""
    output:
        masked = OUTDIR / "masked" / "{sample}.masked.fasta"
    params:
        "pass"
    conda: "../envs/mask.yaml"
    log: "pass"
    shell:"""
    pass
    """