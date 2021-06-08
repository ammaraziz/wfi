suppressMessages(library(dplyr))
suppressMessages(library(stringr))

get_table_data <- function(file_location, org) {
  # returns a list of dataframes read from a list of file names
  if(org == "RSV") {
    dat_out <- list()
    for (f in file_location) {
      df <- read.delim(f, sep = "\t", stringsAsFactors = FALSE)
      if (nrow(df) != 0) {
        splits <- unlist(str_split(string = f, pattern = "/", simplify = F))
        n = length(splits)
        fsplit <- gsub(splits[n], pattern = "-coverage", replacement = "")
        df$sampleID <- paste0(splits[n-2])
        dat_out <- append(dat_out, list(df))
        }
      }
  } else if(org == "FLU") {
    dat_out <- list()
    for (f in file_location) {
      df <- read.delim(f, sep = "\t", stringsAsFactors = FALSE)
      if (nrow(df) != 0) {
        splits <- unlist(str_split(string = f, pattern = "/", simplify = F))
        n = length(splits)
        fsplit <- gsub(splits[n], pattern = "-coverage.txt", replacement = "")
        df$sampleID <- paste0(splits[n-2], '/', fsplit)
        dat_out <- append(dat_out, list(df))
      }
    }
  }
  return(dat_out)
}

calc_stats <- function(coverage_df, org) {
  # coverage_df      : the aa df read in
  # org     : FLU or RSV
  # returns : list(total_mean, total_median, amplicon1_mean, amplicon1_median)
  
  if (toupper(org) == "RSV") {
    out <- coverage_df %>%
      summarise(total_mean = round(mean(Coverage.Depth)), 
                total_median = round(median(Coverage.Depth)),
                amplicon1_mean = round(mean(Coverage.Depth[rsv_primer$loc_start[1]:rsv_primer$loc_end[1]])),
                amplicon1_median = round(median(Coverage.Depth[rsv_primer$loc_start[1]:rsv_primer$loc_end[1]])),
                amplicon2_mean = round(mean(Coverage.Depth[rsv_primer$loc_start[2]:rsv_primer$loc_end[2]])),
                amplicon2_median = round(median(Coverage.Depth[rsv_primer$loc_start[2]:rsv_primer$loc_end[2]])),
                amplicon3_mean = round(mean(Coverage.Depth[rsv_primer$loc_start[3]:rsv_primer$loc_end[3]])),
                amplicon3_median = round(median(Coverage.Depth[rsv_primer$loc_start[3]:rsv_primer$loc_end[3]])),
                amplicon4_mean = round(mean(Coverage.Depth[rsv_primer$loc_start[4]:length(Coverage.Depth)])),
                amplicon4_median = round(median(Coverage.Depth[rsv_primer$loc_start[4]:length(Coverage.Depth)]))
                )
    out$sampleID = coverage_df$sampleID[1]
    out$Manual_Check_Required = abs(out$total_mean/out$total_median - 1) > 0.2
    out$MeanDivMed = abs(out$total_mean/out$total_median - 1)
    out = select(out, sampleID, Manual_Check_Required, MeanDivMed, 1:10)
  }
  
  if (toupper(org) == "FLU") {
      out <- coverage_df %>%
        summarise(total_mean = round(mean(Coverage.Depth)), 
                  total_median = round(median(Coverage.Depth)))
      out$sampleID = coverage_df$sampleID[1]
      out$Manual_Check_Required = abs(out$total_mean/out$total_median - 1) > 0.2
      out$MeanDivMed = abs(out$total_mean/out$total_median - 1)
      out = select(out, sampleID, Manual_Check_Required, MeanDivMed, 1:2)
  }
  
  return(out)
}



#usage
# base_location = "~/Desktop/data/iseq28/"
# 
# locs = get_data_location(base_location, 'cov')
# coverage_data = get_table_data(locs, org = 'FLU')
# 
# out = NULL
# for (i in coverage_data) {
#   tmp = calc_stats(i, org = 'FLU')
#   out = rbind(tmp, out)
# }
# write.table(out, file = "tmp.txt", sep = "\t", row.names = F, quote = F)
