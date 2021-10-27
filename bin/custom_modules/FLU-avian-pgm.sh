# HEADER
PARAM_FILE_NAME="FLU-avian"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-01-28"

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.25				# threshold for insertion refinement
DEL_T=0.60				# threshold for deletion refinement
QUAL_THRESHOLD=19		#average or median threshold for QUALITY reads
MIN_LEN=50
GRID_ON=0

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="LABEL BLAT"
ALIGN_PROG="SAM BLAT"
ASSEM_PROG="SSW"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=19			# average or median threshold for QUALITY reads
MIN_LEN=50				# minimum read length for QUALITY reads
