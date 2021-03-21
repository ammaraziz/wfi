LABEL v0.3.0 by Samuel S. Shepard (vfn4@cdc.gov)
Last updated 2012 May 8, Centers for Disease Control & Prevention.


INPUT
-----
LABEL accepts files in fasta format, whether multi-line or single line and for any number of sequences. Sequences with redundant headers are removed on a first-come first-serve basis. Space characters are also replaced with the underscore character. 


OUTPUT
------
LABEL produces the following output files in a zipped archive: 
	PROJ_final.tab			Standard. Tab-delimited headers & predicted clades.
	PROJ_final.txt			Standard. A prettier output of the above.
	PROJ_LEVEL_trace.tab		Standard. Table of HMM scores at each level, suitable for visualization in R.
	PROJ_results.tab		Standard. For the current LEVEL, tab-delimited headers & predicted clades.
	PROJ_results.txt		Standard. For the current LEVEL, A prettier output of the above.
	FASTA/				Standard. Folder containing fasta files and newick trees.
	FASTA/PROJ_predictions.fas	Standard.  Query sequence file with predictions added like: "_{PRED:CLAD}"
	FASTA/PROT_control.fasta	Optional. Alignment of predictions fasta file and guide sequences.
	FASTA/PROT_control.nwk		Optional. Maximum likelihood tree of the above.
	FASTA/PROJ_reannotated.fas	Default.  Query sequences file with annotations replaced with predicted ones, ordered by clade.
	FASTA/PROJ_ordered.fasta	Optional. Aligned version of the above, still ordered by clade.
	FASTA/PROJ_tree.nwk		Optional. Maximum likelihood tree of the above.
	FASTA/PROJ_clade_CLAD.fas	Standard. The re-annotated file partitioned into separate clade files.
	c-*/				Standard. Clade/lineage subfolder for the hierarchical predictions.
	Read_Me.txt			This file!

*The project name is denoted "PROJ", the prediction level as "LEVEL", the lineage or clade is called "CLAD", and the protein of interest as "PROT".
	

TIME
----
In general, expect no more than 1 second per sequence.  Choosing alignment options may increase the runtime significantly.  Guide sequence libraries are never more than 200 sequences in size.  For the best results using the alignment options, break down the query sequence file into smaller files.


HOW IT WORKS
------------
LABEL or "Lineage Assignment By Extended Learning" uses hidden Markov model (HMM) profiles of lineages/clades--or groups of clades--to score every query sequence and then classify them via machine learning techniques. The HMM scoring step is performed via SAM v3.5 (see http://compbio.soe.ucsc.edu/sam.html for more information).  Scoring is hierarchical in the sense that prediction starts out general (groups of many clades perhaps) and goes to a very specific or terminal level.  This corresponds to the way phylogeny is usually structured.  At each prediction level, data matrices and result files are produced.  The prediction step is done via support vector machines (SVM) using the free SHOGUN Machine Learning Toolbox v1.1.0 (www.shogun-toolbox.org).  Specifically, we utilize a multi-class SVM method called GMNP with an inhomogeneous, normalized polynomial kernel of degree 20.  At the end of the computation, a final output file is produced at the root level with every query sequence being conferred its final prediction.  The input FASTA sequence is annotated with predictions--using the form "{PRED:clade}"--while a second file accepting the predicted annotations as fact is also created.  This second file is sorted by clade and then partitioned into separate files for each clade.  One can optionally align (MUSCLE v3.8.31, see http://www.drive5.com/muscle) the re-annotated FASTA file while retaining clade ordering and then produce a maximum-likelihood tree.  Finally, a guide sequence library containing reference annotations may also be optionally aligned with the query sequences and used to produce a second maximum-likelihood tree for positive control (smaller guide trees may not always be as accurate as fuller ones, when in doubt add more sequences to the clades of interest).  All maximum-likelihood trees use a GTR+GAMMA model with 1000 local support bootstraps and are computed by FastTreeMP v2.1.4 (see http://www.microbesonline.org/fasttree).
