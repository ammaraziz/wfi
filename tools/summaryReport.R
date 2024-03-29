#!/usr/bin/env Rscript

# this script runs from the command line
# it generates a report from an IRMA run to summaries all the genes/samples


library(optparse, quietly=TRUE)

option_list <- list(
  make_option(c("-i","--input"), help = "Input base directory", 
              action="store", type="character", default=NA),
  make_option(c("-o","--output"),help = "Output directory", 
              action="store", type="character", default=NA),
  make_option(c("-r","--organism"), help = "Name of organism (flu or rsv)", 
              action="store", type="character", default=NA),
  make_option(c("-n","--verbose"),help = "Be verbose", 
              action="store_true", default=FALSE)
)
parser <- OptionParser(
  usage = paste("%prog -i [INPUTDIR] -o [OUTPUTDIR] -r [FLU or RSV]",
                "This script will parse logs of IRMA to create an easy to read pdf", sep="\n"),
  epilogue = "INPUT, OUTPUT and ORG are required",
  option_list=option_list)

#custom function to stop quietly
stop_quietly = function(message) {
  opt = options(show.error.messages = FALSE)
  on.exit(options(opt))
  cat(message, sep = "\n")
  quit()
}

usage_message = "Failed to parse command-line parameters.

This script parses IRMA output files in a directory to produce useful plots collated in one location. 
Usage:
      summaryReport.R -i [INPUTDIR] -o [OUTPUTDIR] -r [FLU or RSV]
      summaryReport.R -i .../assemblies/ -o .../report/ -r flu

Input directory needs to be the wfi output directory called 'assemblies'
The -r or --organism flag is needed because IRMA produces a file for each assembled contig (for flu thats 8 genes)
"

cov_message = "# Manual_Check_Required column will indicate (TRUE) if the depth is more than 20% different between the mean and median calculation. [mean(reads)/median(reads) > 20%]"
qc_message = "# qc_pass_50pc  - Did 50% of the total reads pass IRMA QC? [total_reads_pass_qc > 50%]
# nomatch_high - Num of no-match reads that is greater than match-reads. [nomatch > match]
# match_adequate - Num of match-reads above a minimum threshold (approx. ~50x depth). [reads_match >= min_read]
# alt_metric - Large Num of alt-reads detected - helps to identify mixture (10% of reads and approx. 25x depth). [altmatch/match > 0.10 & altmatch > min_read_alt] "

arguments=NA
tryCatch(
  { arguments = parse_args(object = parser, positional_arguments = TRUE) },
  error = function(e) { })
if (any(is.na(arguments$options))) {
  stop_quietly(message = usage_message)
  print(parser$usage)
}
opts = arguments$options

## Parameter Validation
if ( is.na(opts$input) ) { stop("Missing --input flag") }
if ( is.na(opts$output) ) { stop("Missing --output flag") }
if ( is.na(opts$organism) ) { stop("Missing --organism flag. Options: RSV or FLU") }
if ( !(toupper(opts$organism) %in% c('RSV', 'FLU')) ) {
  print(opts)
  stop(paste("Organism selected:", opts$organism, 'is not RSV or FLU.'))}

###################### vars for plotting
rsv_gene_locs = read.table("tools/colors/rsv.tsv", header = T, sep = "\t")
rsv_primer = data.frame(start = c(49, 3944, 7215, 10959),
                        end = c(4049, 7528, 11165, 15333))
flu = read.table("tools/colors/flu.tsv", header = T, sep = "\t", comment.char = "")
gcolor = split(flu$colors, flu$gene)
# gcolor = list("A_HA_H1" = '#a6cee3', 
#               "A_HA_H3" = '#a6cee3', 
#               "A_HA_H5" = '#a6cee3',
#               "A_HA_H7" ='#a6cee3',
#               "A_HA_H9" ='#a6cee3',
#               "A_HA_H10" = '#a6cee3',
#               "A_MP" = '#1f78b4',
#               "A_NA_N1" = '#b2df8a',
#               "A_NA_N2" = '#b2df8a',
#               "A_NA_N4" = '#b2df8a',
#               "A_NA_N5" = '#b2df8a',
#               "A_NA_N6" = '#b2df8a',
#               "A_NA_N7" = '#b2df8a',
#               "A_NA_N8" = '#b2df8a',
#               "A_NP" = '#33a02c',
#               "A_NS" = '#fb9a99',
#               "A_PA" = '#ff7f00',
#               "A_PB1" = '#6a3d9a',
#               "A_PB2" = '#b15928',
#               "B_HA" = '#a6cee3',
#               "B_NA" = '#b2df8a',
#               "B_NP" = '#33a02c',
#               "B_NS" = '#fb9a99',
#               "B_PA" = '#ff7f00',
#               "B_PB1" = '#6a3d9a',
#               "B_PB2" = '#b15928',
#               "B_MP" = '#1f78b4')



###################################

source("./tools/depthPlots.R")
source("./tools/depthTable.R")
source("./tools/qcTable.R")

main  <- function() {

  opts$organism = toupper(opts$organism)
  
  # get file names
  base_location = as.character(opts$input)
  file_names = file_finder(base_location)
  
  # calc average depth
  data_aa = get_data(file_location = file_names$alleleDepth, file_name = file_names$full_name)
  data_aa_df = data.frame(full_name = names(data_aa))
  data_aa_df$data_aa = data_aa
  file_names = left_join(file_names, data_aa_df, by = "full_name")
  # clean up
  rm(data_aa, data_aa_df)
  
  # qc table
  qc_loc = get_data_location(base_location = base_location, type = 'rc')
  qc_data = get_qc_data(qc_loc)
  qc_stats = calc_qc_stats(qc_data, org = opts$organism)
  write.table(qc_stats, file = paste0(opts$output, 'summary_QC_table.tsv'), sep = "\t", row.names = F, quote = F)
  write(qc_message, file = paste0(opts$output, 'summary_QC_table.tsv'), append = T)
  
  # coverage table
  locs = get_data_location(base_location = base_location, type = 'cov')
  coverage_data = get_table_data(locs, org = opts$organism)
  out_table = NULL
  for (i in coverage_data) {
    tmp = calc_stats(i, org = opts$organism)
    out_table = rbind(tmp, out_table)
  }
  write.table(out_table, file = paste0(opts$output, 'summary_depth_table.tsv'), sep = "\t", row.names = F, quote = F)
  write(cov_message, file = paste0(opts$output, 'summary_depth_table.tsv'), append = T)
  
  if (opts$organism == "RSV") {# RSV
    final_out = plot_combine_rsv(file_names)
  } else if (opts$organism == "FLU") {# FLU
    final_out = plot_combine_flu(file_names)
  } else {
    stop("Error: Organism not RSV or FLU. Check input")
  }
  ggsave(filename =  paste0(opts$output, "depth_coverage_report.pdf"),
         plot = final_out, device = "pdf",
         dpi = 300,  width = 420, height = 297, units = 'mm')
}

main()

quit()
