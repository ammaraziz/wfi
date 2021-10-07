# HEADER
PARAM_FILE_NAME="fast"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-02-27"

# READ GATHERING & ASSEMBLY
NO_MERGE=1		# no read pair merging
MAX_ROUNDS=1		# max read gather rounds
MAX_ITER_SSW=1		# max SSW iterations to perform

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="BLAT"
ASSEM_PROG="SSW"
