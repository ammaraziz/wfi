rule makeBedForMasking:
    input:
        irmastatus = rules.irma.output.status,
        irmadir = OUTDIR / "irma" / "{sample}/",
    output:
        bed = OUTDIR / "mask" / "{sample}" / "{sample}.mask.bed",
        status = OUTDIR / "status" / "mask_{sample}.txt",
    params:
        mincov = 20
    run:
        irma_files = get_irma_files(input['irmadir'])
        make_bed_for_masking(infile=irma_files['coverage'][0],
                             chromCol=0,
                             depthCol=6,
                             outfile=output['bed'],
                             min_cov=params.mincov)
        with open(output['status'], "w") as file:
            file.write("")

# rule maskFasta:
#     input:
#         bed = rules.makeBedForMasking.output.bed
#     output:
#         maskedfasta = OUTDIR / "mask" / "{sample}" / "{sample}.masked.fasta"
#     threads: 1
#     conda: ""
#     params:
#         IRMA = OUTDIR
#     shell:"""
#     bedtools maskfasta \
#     fi {input.fasta} \
#     -bed {input.bed} \
#     -fo {output.maskedfasta}
#     """