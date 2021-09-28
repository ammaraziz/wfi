suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

min_reads_metric = c('RSV' = 5000, 'FLU' = 2500)

get_qc_data <- function(file_location) {
  # returns a list of dataframes read from a list of file names
  dat_out <- data.frame()
  for (f in file_location) {
    df = read.delim(f, sep = "\t", stringsAsFactors = FALSE)
    df = df %>% select(Record, Reads) 
    if (nrow(df) != 0) {
      splits <- unlist(str_split(string = f, pattern = "/", simplify = F))
      n = length(splits)
      fsplit <- gsub(splits[n], pattern = "-coverage", replacement = "")
      df$sampleID <- paste0(splits[n-2])
      df_wide = pivot_wider(df, names_from = Record, values_from = Reads)  
      dat_out <- bind_rows(dat_out, df_wide)
    }
  }
  names(dat_out) = make.unique(gsub("\\d-", "", names(dat_out)))
  return(dat_out)
}

calc_qc_stats = function(dat_out, org) {
  # dat_out   :   wide dataframe output from get_qc_data
  # returns dat_out with additional QC columns
  min_read = min_reads_metric[[org]] * 2
  min_read_alt = min_reads_metric[[org]]
  
  qc_pass_50pc = (dat_out$passQC/dat_out$initial) > 0.50
  nomatch_high = (dat_out$match/dat_out$passQC) > (dat_out$nomatch/dat_out$passQC)
  match_adequate = (dat_out$match + dat_out$altmatch) >= min_read
  alt_metric = ((dat_out$altmatch/dat_out$match) > 0.10) & (dat_out$altmatch > min_read_alt)
  
  dat_out = dat_out %>% mutate(qc_pass_50pc = case_when(qc_pass_50pc == TRUE ~ 'ok',
                                                        qc_pass_50pc == FALSE ~ 'FAILED',
                                                        TRUE ~ NA_character_),
                               nomatch_high = case_when(nomatch_high == TRUE ~ 'FAILED',
                                                        nomatch_high == FALSE ~ 'ok',
                                                        TRUE ~ NA_character_
                                                        ),
                               match_adequate = case_when(match_adequate == TRUE ~ 'ok',
                                                          match_adequate == FALSE ~ 'FAILED',
                                                          TRUE ~ NA_character_),
                               alt_metric = case_when(alt_metric == TRUE ~ 'Mix?',
                                                      alt_metric == FALSE ~ '-',
                                                      TRUE ~ NA_character_)
  )
  dat_out = dplyr::select(dat_out, sampleID, qc_pass_50pc, nomatch_high, match_adequate, alt_metric, 2:(ncol(dat_out) - 4))
  return(dat_out)
}


# # usage
# file = ".../tables/READ_COUNTS.txt"
# 
# dat = read.delim(file, sep = "\t")
# a = get_qc_data(file)
# b = calc_qc_stats(a, org = 'RSV')
 