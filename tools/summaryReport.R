#!/usr/bin/env Rscript

# this script runs from the command line
# it generates a report from an IRMA run to summaries all the genes/samples


library(optparse,quietly=TRUE)

option_list <- list(
  make_option(c("-i","--input"), help = "Input base directory", 
              action="store", type="character", default=NA),
  make_option(c("-o","--output"),help = "Output directory", 
              action="store", type="character", default=NA),
  make_option(c("-r","--organism"), help = "Name of organism (flu or rsv)", 
              action="store", type="character", default=NA),
  make_option(c("-v","--verbose"),help = "Be verbose", 
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
rsv_gene_locs = list(ns1 = 576, ns2 = 1097, n = 2329, p = 3220, m = 4180, sh = 4599, g = 5565, f = 7521, m2 = 8527, l = 15037)

rsv_primer = list(loc_start = c(49, 3944, 7215, 10959),
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
              "A_PB2" = '#b15928')


###################################

source("./tools/depthPlots.R")

main  <- function() {

  #print params for user
  write(paste0("organism: ", opts$organism), stdout())
  write(paste0("input: ", opts$input), stdout())
  write(paste0("output: ", opts$output), stdout())

  base_location = as.character(opts$input)

  locations_aa = get_data_location(base_location, type = 'aa')
  data_aa = get_data(locations_aa)
  data_aa_avg = lapply(data_aa, FUN = average_counts, type = 'aa', by_num = 5)

  #RSV
  if (opts$organism == "RSV") {
    names(data_aa_avg) = str_match(locations_aa, pattern = "(\\w+)/tables")[,2]
    #plot
    final_out = plot_combine_rsv(data_avg = data_aa_avg, base_location = base_location)
  
    #FLU
  } else if (opts$organism == "FLU") {
    # rename
    file_names = tibble(
      sample = str_match(locations_aa, pattern = "(N100.+?)/tables")[, 2],
      file_name = str_match(locations_aa, pattern = "[A|B]_\\w{2,3}.+")) %>% 
      mutate(gene = str_replace(file_name, pattern = "-allAlleles.txt", replacement = "")) %>%
      mutate(sample_gene = paste0(sample, "/", gene))
    
    names(data_aa_avg) = file_names$sample_gene
    #plot
    final_out = plot_combine_flu(data_avg = data_aa_avg,
                             gene = file_names$gene,
                             sample = file_names$sample,
                             base_location = base_location)
  } else {
    stop("Error: Organism not RSV or FLU. Check input")
  }
  ggsave(filename =  opts$output,
         plot = final_out, device = "pdf",
         dpi = 300,  width = 420, height = 297, units = 'mm')
}

main()

quit()
