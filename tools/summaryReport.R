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
      tmp_script.R -i [INPUTDIR] -o [OUTPUTDIR] -r [FLU or RSV]
      tmp_script.R -i .../assemblies/ -o .../report/ -r flu

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
rsv_gene_locs = data.frame(loc_start = c(44, 596, 1125, 2331, 3224, 4190, 4644, 5619, 7567, 8460),
                           loc_end = c(576, 1097, 2329, 3220, 4180, 4599, 5565,7521, 8527, 15037),
                           gene_name = c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"))

rsv_primer = data.frame(loc_start = c(49, 3944, 7215, 10959),
                      loc_end = c(4049, 7528, 11165, 15333))

gcolor = list("A_HA_H1" = '#a6cee3',
              "A_HA_H3" = '#a6cee3',
              "A_MP" = '#1f78b4',
              "A_NA_N1" = '#b2df8a',
              "A_NA_N2" = '#b2df8a',
              "A_NP" = '#33a02c',
              "A_NS" = '#fb9a99',
              "A_PA" = '#ff7f00',
              "A_PB1" = '#6a3d9a',
              "A_PB2" = '#b15928',
              "B_HA" = '#a6cee3',
              "B_N2" = '#b2df8a',
              "B_NP" = '#33a02c',
              "B_NS" = '#fb9a99',
              "B_PA" = '#ff7f00',
              "B_PB1" = '#6a3d9a',
              "B_PB2" = '#b15928',
              "A_MP" = '#1f78b4')



###################################

source("./tools/depthPlots.R")
source("./tools/depthTable.R")
source("./tools/qcTable.R")

main  <- function() {

  #print params for user
  write(paste0("organism: ", opts$organism), stdout())
  write(paste0("input: ", opts$input), stdout())
  write(paste0("output: ", opts$output), stdout())

  base_location = as.character(opts$input)

  locations_aa = get_data_location(base_location, type = 'aa')
  data_aa = get_data(locations_aa)
  data_aa_avg = lapply(data_aa, FUN = average_counts, type = 'aa', by_num = 5)

  # qc table
  qc_loc = get_data_location(base_location = base_location, type = 'rc')
  qc_data = get_qc_data(qc_loc)
  qc_stats = calc_qc_stats(qc_data, org = 'RSV')#opts$organism)
  write.table(qc_stats, file = paste0(opts$output, '_qcTable.tsv'), sep = "\t", row.names = F, quote = F)
  write(qc_message, file = paste0(opts$output, '_qcTable.tsv'), append = T)
  
  # RSV
  if (opts$organism == "RSV") {
    names(data_aa_avg) = str_match(locations_aa, pattern = "(\\w+)/tables")[,2]
    # plot
    final_out = plot_combine_rsv(data_avg = data_aa_avg, base_location = base_location)
    # coverage table
    locs = get_data_location(base_location = base_location, type = 'cov')
    coverage_data = get_table_data(locs, org = 'RSV')
    out_table = NULL
    for (i in coverage_data) {
      tmp = calc_stats(i, org = 'RSV')
      out_table = rbind(tmp, out_table)
    }
    write.table(out_table, file = paste0(opts$output, '_depthTable.tsv'), sep = "\t", row.names = F, quote = F)
    write(cov_message, file = paste0(opts$output, '_depthTable.tsv'), append = T)
  
    # FLU
  } else if (opts$organism == "FLU") {
    # rename
    file_names = tibble(
      #sample = str_match(locations_aa, pattern = "(\\w+)/tables")[, 2],
      sample = str_match(locations_aa, pattern = "([A-Za-z0-9_-]+_S\\d{1,3})\\/(?!\\/tables)")[, 2],
      file_name = str_match(locations_aa, pattern = "[A|B]_\\w{2,3}.+")) %>% 
      mutate(gene = str_replace(file_name, pattern = "-allAlleles.txt", replacement = "")) %>%
      mutate(sample_gene = paste0(sample, "/", gene))
    
    names(data_aa_avg) = file_names$sample_gene
    ##plot
    final_out = plot_combine_flu(data_avg = data_aa_avg,
                             gene = file_names$gene,
                             sample = file_names$sample,
                             base_location = base_location)
    ## cov table
    locs = get_data_location(base_location = base_location, type = 'cov')
    coverage_data = get_table_data(locs, org = 'FLU')
    out_table = NULL
    for (i in coverage_data) {
      tmp = calc_stats(i, org = 'FLU')
      out_table = rbind(tmp, out_table)
    }
    write.table(out_table, file = paste0(opts$output, '_depthTable.tsv'), sep = "\t", row.names = F, quote = F)
    write(cov_message, file = paste0(opts$output, '_depthTable.tsv'), append = T)
  } else {
    stop("Error: Organism not RSV or FLU. Check input")
  }
  ggsave(filename =  paste0(opts$output, "depthReport.pdf"),
         plot = final_out, device = "pdf",
         dpi = 300,  width = 420, height = 297, units = 'mm')
}

main()

quit()
