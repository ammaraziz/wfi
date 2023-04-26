## wfi (WhoFlu IRMA) 

Snakemake pipeline for running IRMA. Designed for Influenza and RSV Illumina Sequencing.

Note: This pipeline is not ready for general use, while it generally works, it is very brittle.

### Requirements:
- Linux Distro (or \*unix system like MacOS) 
- Conda
- Snakemake
- Cutadapt
- IRMA
- R
- R Packages: `ggplot2` `dplyr` `stringr` `tidyr` `cowplot` `gridExtra` `furrr`

---

### Installation

1. Install `miniconda`: https://docs.conda.io/en/latest/miniconda.html
2. Install `snakemake`: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
    ```
    conda install -n base -c conda-forge mamba
    mamba install -c bioconda -c conda-forge snakemake-minimal python=3.9
    ```
3. Install `cutadapt` and `biopython`:
    ```
    mamba install -c bioconda -c conda-forge cutadapt biopython bbmap
    ```
4. Install R (>3.6 should work) and R packages:
    ```
    mamba install -c conda-forge r-base 
    mamba install -c r r-ggplot2 r-dplyr r-tidyr r-cowplot r-gridExtra r-optparse r-furrr
    ```
 Note: installing r packages through conda is troublesome for some, if so install manually in R. 
   
5. Install custom verison of IRMA which contains the RSV module:
    ```
    mamba install -c ammaraziz irma 
    ```
6. Finally, download this repository and store in your /bin/


---
### Usage

To use the pipeline, follow these steps:

1. Navigate to `wfi_config.yaml` and modify as appropriate:

| Params | Values | Information |
| --- | --- | --- |
| input_dir | path | input directory - location of the raw fastq files for input |
|output_dir|path|output directory - location to output results - same dir where the config sits|
|second_assembly|`True`/`False`|if you suspect mixtures, set to `True` . It will increase run time substantially|
|subset|`True`/`False`|if you are only sequencing HA/NA/MP set this to `True` else leave as `False`|
|trim_prog|`standard`/`tile`| Trimming program to use, tile (bbduk) or standard (cutadapt)|
|trim_org|`h1`/`h3`| Influenza only, Flu subtype|
|technology|`illumina`/`ont`/`pgm`| seq technology used, will change the module by IRMA|


2. Check snakemake is installed, if an error is produced it means snakemake was not found or it is not installed.
```	
% snakemake --version
% 5.10.0 
```

3. Test the pipeline, this will output all the commands that will be run. Look for errors (red).
```
% snakemake -np
```

4. Run the pipeline, with option -j to specify number of cores to use.
```
% snakemake -j 8
```
---
### Output structure:
	
1. Pipeline will output correctly formatted names located in:
	
> `{output_dir}`/assemblies/rename/

2. Sorted by subtype - most likely the disired output:
	
> `{output_dir}`/assemblies/rename/type/FLU`{A|B}`

3. IRMA assembly specific files, see: https://wonder.cdc.gov/amd/flu/irma/output.html
	
> `{output_dir}`/assemblies/`{sampleID}`/

4. Files for depth and summary info located in:
	
> `{output_dir}`/assemblies/`{sampleID}`/figures/
> `{output_dir}`/assemblies/`{sampleID}`/tables/
	


## Dependencies

- [BLAT] for the match step
- [LABEL], which also packages certain resources used by IRMA:
    - Sequence Alignment and Modeling System ([SAM]) for both the rough align and sort steps
    - [Shogun Toolbox], which is an essential part of [LABEL], is used in the sort step
- [SSW] for the final assembly step, download our [minor modifications to SSW](https://wonder.cdc.gov/amd/flu/irma/ssw-modified.tar.gz)
- [samtools] for BAM-SAM conversion as well as BAM sorting and indexing
- [GNU Parallel] for single node parallelization
- R and these R packages: optparse, ggplot2, dplyr, tidyr, stringr, cowplot, gridExtra
- 
[IRMA manuscript]: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6
[IRMA Website]: https://wonder.cdc.gov/amd/flu/irma/
[IRMA v0.6.7 zip file]: https://wonder.cdc.gov/amd/flu/irma/flu-amd-201704.zip
[SAM]: https://www.ncbi.nlm.nih.gov/pubmed/9927713
[BLAT]: http://www.kentinformatics.com/products.html
[LABEL]: https://wonder.cdc.gov/amd/flu/label/
[SSW]: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138
[samtools]: http://www.htslib.org/
[GNU Parallel]: https://www.gnu.org/software/parallel/
[Shogun Toolbox]: http://shogun.ml/	
	
---
### Troubleshooting problems:
```
1. Error regarding path directories		Check input and output directorys you've specified end with a '/'

2. Error: Nothing to be done			Check config file and ensure you've changed the input/output directories. 

3. A job crashed. What do I do?			Two options, delete the output directory so snakemake can run everything again. 
						Or find out where it crashed and delete the whole folder/sample. 
						Example, sometimes IRMA produces errors, find the sample which crashed, 
						go to assemblies and delete the corresponding {sampleID} folder. Rerun snakemake.
4. I'm very confusd 										
   or I need more help			
   or I've screwed something up badly!		Shoot me an email
```

For any issues please submit a github issue.
