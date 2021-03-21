WFI (WhoFlu IRMA) is a pipeline that uses Snakemake to handle the job submission. It is designed only for Influenza Illumina Sequencing.

To use the pipeline, follow these steps:

	1. Copy wfi_config.yaml to a new directory
	2. Modify wfi_config.yaml:
		~ input_dir 	 		input directory - location of the raw fastq files for input.
		~ output_dir 	 		output directory - location to output results - same dir where the config sits.
		~ subset 		 		if you are only sequencing HA/NA/MP set this to True else leave as False.
		~ second_assembly		if you suspect mixtures, set to True . It will increase run time substantially.

	3. Check snakemake is installed, if an error is produced it means snakemake was not found or it is not installed.
		%	snakemake --version
		output:	5.10.0

	4. Test the pipeline, this will output all the commands that will be run. Look for errors (usually in red).
		% snakemake -np

	5. Run the pipeline, with option -j to specify number of cores to use (recommended is no more than 46)
		% snakemake -j 46


Output structure:
	1. Pipeline will output correctly formatted names located in:
		{output_dir}/assemblies/rename/

	2. Sorted by subtype - most likely the disired output:
		{output_dir}/assemblies/rename/type/FLU{A|B}

	3. IRMA assembly specific files, see: https://wonder.cdc.gov/amd/flu/irma/output.html
		{output_dir}/assemblies/{sampleID}/

	4. Files for depth and summary info located in:
		{output_dir}/assemblies/{sampleID}/figures/
		{output_dir}/assemblies/{sampleID}/tables/

Troubleshooting problems:

	1. Error regarding path directories		Check input and output directorys you've specified end with a '/'

	2. Error: Nothing to be done			Check config file and ensure you've changed the input/output directories.

	3. A job crashed. What do I do?			Two options, delete the output directory so snakemake can run everything again.
											Or find out where it crashed and delete the whole folder/sample.
											Example, sometimes IRMA produces errors, find the sample which crashed,
											go to assemblies and delete the corresponding {sampleID} folder. Rerun snakemake.

	4. I'm very confusd
	or I need more help						Contact Ammar via email: ammar.aziz@mh.org.au or go bug him in person for help at any time.
	or I've screwed something up badly!