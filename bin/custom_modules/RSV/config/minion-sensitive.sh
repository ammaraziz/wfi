# HEADER
PARAM_FILE_NAME="CoV MinION Sensitive"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.1"
PARAM_FILE_DATE="2021-02"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=0	# average or median threshold for QUALITY reads
MIN_LEN=125		# minimum read length for QUALITY reads
INS_T=0.75		# threshold for insertion refinement
DEL_T=0.90		# threshold for deletion refinement : 1 => turn OFF deletion editing
MIN_RP=1		# minimum read pattern count to continue
MIN_RC=3		# minimum read count to continue
MIN_CONS_SUPPORT=3

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8		# minimum average variant quality, does not apply to deletions
MIN_FI=0.01		# minimum insertion variant frequency
MIN_FD=0.02		# minimum deletion variant frequency

ALIGN_PROG="SAM"	# rough alignment with HMM
ASSEM_PROG="SSW"	# final assembly with SSW

SSW_M=2			# smith-waterman match score
SSW_X=3			# smith-waterman mismatch penalty
SSW_O=6			# smith-waterman gap open penalty
SSW_E=1			# smith-waterman gap extension penalty
