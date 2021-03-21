          SAM:  Sequence Alignment and Modeling Software System
                              Binary Distribution

This is the README for the binary distribution of SAM.  This is copyrighted
software, please see the file COPYRIGHT in this directory for details.
For more information about SAM, visit the SAM WWW site at:

        http://www.soe.ucsc.edu/research/compbio/sam.html

Distribution organization:

    o arch/bin/* - SAM programs and scripts for the system
      architecture.  arch is a canonical system name that contains the
      hardware architecture and operating system version of the system
      that SAM was built on.  There normally does not have to be an
      exact operating system version match to run the programs,
      however this information maybe useful if compatibility problems
      are encountered.  See below for further discussion of SAM-T99
      installation.

    o arch/lib/sam/* - SAM runtime files.  This includes the prior
      libraries, reguarlizers, colors files and (as a subdirectory)
      Smith and Waterman cost matrices.
      These files are normally installed in /usr/local/lib/sam/*.  
      ***** If installed elsewhere, the environment variable
            PRIOR_PATH must be set to that directory.
      ***** If the matrix directory is not installed as a subdirectory
            of PRIOR_PATH, set the BLASTMAT environment variable to
	    the path of the matrix directory for performing SAM Smith
	    and Waterman.

    o arch/doc/sam/* - HTML, PDF,  and postscript documentation.

    o sam-demos - Demonstration programs.
      ***** Modify the path names at the top of the Makefile in the
	    sam-demos directory to match your setup.
      ***** Be sure to use complete (non-relative) path names

The sae and hmmedit programs are in a separate distribution (samtools)
that contains an `arch' directory as described above as well as a
patch that is required for rasmol to work with sae.


Installation:
    o General installation

    Move the files to the desired directories on your system.

    o SAM/Target2k configuration	    

    In the final bin directory, run:

    ./sam-t-configure `which perl`

    to set the perl path for PERL5 (depending on your local setup,
    you may need to use `which perl5` or type in the pathname
    directly), and then edit sam-t2k.conf (and sam-t99.conf) as
    follows:  

    ncbi_blast_prog:    full path of the NCBI blastall program
			change to the data directory for psiblast.
			With appropriate psiblast installation, you
			may not	need this.

    wu_blast_prog:      If you wish to use wu-blastp, put its path
			name here, and use -wublast on the target2k
			command line

    NR :		A large non-redunundant database is required
			for SAM-T99, and must be blast-indexed.  We
			use the nr from NCBI (200MB Dec 2001):
			ftp://ncbi.nlm.nih.gov/blast/db/nr.Z 

    sam_bin_dir
    sam_t99_bin_dir:   Set to where you have installed the SAM
		       programs and scripts.
    sam_t2k_bin_dir:   Set to where you have installed the SAM
		       programs and scripts.
    reg_lib_dir:       Set to where you have installed the prior
		       libraries. 

    gunzip_path:       Make sure that the GNU gunzip programs are in
		       this path.

    unix_path:	       Where standard unix programs are.		       


    Also note that target99 and target2k will automatically give full
    group access to all of its created files and directories, and will
    change the group to be 'protein' if such a group exists in the
    user's group list as a result of executing the "groups" command.
    If you have a different default group, change SamT99.pm and SamT2K.pm
    make_directories routine.
    







