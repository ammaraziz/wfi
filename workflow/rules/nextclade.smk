rule nextclade_sort:
    message: "Running nextclade sort to detect best dataset"
    input: 
        pass
    output:
        results =  WORKDIR / "nextclade" / "result.dataset.tsv"
    params:
        pass
    threads: config['threads']['nextclade']
    log: WORKDIR / "logs" / "nextclade.sort.txt"
    conda: "../envs/qc.yaml"
    shell:"""
    nextclade sort \
    --quiet \
    --jobs {threads} \
    --output-results-tsv {output.dataset_results} \
    {input.consensus} 2> {log}
    """

rule nextclade_split:
    message: "Splitting Nextclade Sort Results"
    input:
        rules.nextclade_sort.output.results
    output:
        pass
    params:
        pass
    threads: 1
    log: WORKDIR / "logs" / "nextclade.split.txt"
    conda: "../envs/qc.yaml"
    shell:"""
    echo test
    """

rule nextclade_run:
    message: "Running Nextclade on samples"
    input:
    output:
    params:
    threads: 10
    log: WORKDIR / "logs" / "nextclade.run.txt"
    conda: "../envs/qc.yaml"
    shell:"""
    echo test
    """