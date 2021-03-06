### there is an issue with the two plot functions for RSV/FLU.
### They are not harmonised, i.e. they have different inputs. Make sure the inputs are identical for easy use.

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(cowplot))
suppressMessages(library(gridExtra))

options(dplyr.summarise.inform = FALSE)
options(warn = -1)
debuggingState(on = FALSE)


get_data_location <- function(base_location, type) {
  # searches a given directory and returns names table files generated by IRMA.
  pats <- list(
    rc = "READ_COUNTS.txt", del = "*-deletions.txt", var = "*-variants.txt",
    aa = "*-allAlleles", ins = "*-insertions.txt", cov = "*-coverage.txt",
    ps = "*-paringStats.txt"
  )

  if (type %in% pats) {
    message("Error - Could not find tables to parse from pattern: " + paste0(type))
  } else {
    pattern <- pats[[type]]
  }
  f <- list.files(base_location, pattern = pattern, full.names = TRUE, recursive = TRUE)
  return(f)
}


get_minor_vars <- function(file_location, ...) {
  # reads in vcf file to get minor variants
  # return df of: pos ref alt dp af
  extra <- list(...)
  tryCatch(
    vcf <- read.delim(file_location,
      header = F, skip = 21, sep = "\t",
      col.names = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"),
      colClasses = "character"
    ),
    error = function(e) {
      message(paste0("Error, VCF file doesn't exist in", file_location))
      return(0)
    }
    # finally = {message("Error, VCF file doesn't exist");return(0)}
  )
  if (nrow(vcf) == 0) {
    #write("VCF file is empty \n", stdout())
    return(0)
  }

  vcf <- vcf %>%
    filter(filter == "PASS") %>%
    separate(
      col = info, into = c("dp", "af", "aq", "pub", "qub"),
      sep = ";", extra = "merge"
    )
  vcf$dp <- as.numeric(gsub(x = vcf$dp, pattern = "DP=", replacement = ""))
  vcf$af <- as.numeric(gsub(x = vcf$af, pattern = "AF=", replacement = ""))
  vcf$pos <- as.numeric(vcf$pos)

  output <- vcf[, c("pos", "ref", "alt", "dp", "af")]
  return(output)
}

get_data <- function(file_location) {
  # returns a list of dataframes read from a list of file names

  dat_out <- list()
  for (f in file_location) {
    df <- read.delim(f, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(df) != 0) {
      splits <- unlist(str_split(string = f, pattern = "/", simplify = F))
      fsplit <- gsub(splits[8], pattern = "-allAlleles.txt", replacement = "")
      df$sampleID <- paste0(splits[6], "/", fsplit)
      dat_out <- append(dat_out, list(df))
    }
  }
  return(dat_out)
}

average_counts <- function(df, type, by_num) {
  # type    : the aforementioned file patterns: rc, del, var, aa, ins, cov, ps
  # by_num  : the number of rows to average by
  # returns : a dataframe similar to input but averaged over by_num

  ## type = aa
  if (type == "aa") {
    df_out <- df %>%
      filter(Allele_Type == "Consensus") %>%
      group_by(group_var = trunc(2:(n() + 1) / by_num)) %>%
      summarise(coverage = mean(Count), genomic_position = min(Position))
  } else {
    return("this is tmp, probably an error")
  }
  return(df_out)
}

plot_flu_cov <- function(df, sample, gene, base_location) {
  # plots coverage given a dataframe
  # df        :   dataframe  - output of average_counts function
  # plot_name :   character of plot name
  
  vcf <- list.files(
    path = paste(base_location, sample, sep = "/"),
    pattern = paste0(gene, ".vcf"), full.names = T
  )

  minor_vars <- get_minor_vars(vcf)
  plot_name <- paste(sample, gene, sep = "/")

  # catch empty or missing vcf files
  if (minor_vars == 0 | minor_vars == 00) {
    annotations <- list()
  } 
  else {
    minor_loc <- minor_vars$pos
    minor_dp <- round(minor_vars$dp * minor_vars$af, digits = 1)
    snp <- paste0(minor_vars$ref, minor_vars$pos, minor_vars$alt, " :DP", minor_dp)

    y_centre <- mean(range(df$coverage))
    y_max <- max(df$coverage) - ((y_centre / 20) * 3)
    x_minor_unit <- mean(range(df$genomic_position)) / 30

    annotations <- list(
      geom_vline(xintercept = minor_loc, color = "grey"),
      annotate(
        geom = "text", x = minor_loc,
        y = y_max, label = snp, color = "black", angle = 90, size = 3
      )
    )
  }

  p <- ggplot(df, aes(x = genomic_position, y = coverage)) +
    geom_col(fill = gcolor[[gene]]) +
    annotations +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle(paste(plot_name)) +
    coord_cartesian(clip = "off")

  return(p)
}

plot_rsv_cov <- function(df_avg, plot_name, base_location) {
  # plots coverage given a dataframe.
  # df          :   dataframe  - output of average_counts()
  # plot_name   :   character of plot name

  # catch empty or missing vcf files

  vcf <- list.files(
    path = paste(base_location, plot_name, sep = "/"),
    pattern = ".vcf", full.names = T
  )

  minor_vars <- get_minor_vars(vcf)

  if (minor_vars == 0 | minor_vars == 00) {
    annotations <- list()
  } else {
    minor_loc <- minor_vars$pos
    minor_dp <- round(minor_vars$dp * minor_vars$af, digits = 1)
    snp <- paste0(minor_vars$ref, minor_vars$pos, minor_vars$alt, " - DP:", minor_dp)

    y_centre <- mean(range(df_avg$coverage))
    y_max <- max(df_avg$coverage) - ((y_centre / 20) * 3)
    x_minor_unit <- mean(range(df_avg$genomic_position)) / 30

    annotations <- list(
      geom_vline(xintercept = minor_loc, color = "black", size = 0.25),
      annotate(
        geom = "text", x = minor_loc,
        y = y_max, label = snp, color = "black", angle = 90, size = 2.5))
  }

  p <- ggplot(df_avg, aes(x = genomic_position, y = coverage)) +
    # primer regions
    geom_rect(data = rsv_primer, inherit.aes = FALSE,
              aes(xmin = loc_start,
                  xmax = loc_end,
                  ymin = -Inf, ymax = Inf, 
                  fill = c("red", "green", "blue", "yellow")),
              alpha = 0.3) +
    
    geom_col() +
    geom_hline(yintercept = -50) +
    # genes displayed below
    geom_rect(data = rsv_gene_locs, inherit.aes = FALSE,
              aes(xmin = loc_start,
                  xmax = loc_end, 
                  ymin = 0, ymax = -100), 
              color = c("red", "grey", "green", "yellow", "red", "grey", "green", "yellow", "red", "grey"),
              fill = c("red", "grey", "green", "yellow", "red", "grey", "green", "yellow", "red", "grey")) +
    # gene names
    geom_text(data = rsv_gene_locs,
              aes(x = loc_start + (0.05 * loc_start), 
                  y = -50, 
                  label = gene_name), 
              size = 3) +

    annotations +
    
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle(paste(plot_name)) +
    coord_cartesian(clip = "off") + 
    theme(legend.position = "none")

  return(p)
}

plot_combine_rsv <- function(data_avg, base_location) {
  # combines plots into one neat page
  # data_avg    :   named list of averaged and combined list of data

  plots <- list()
  for (i in seq_along(data_avg)) {
    p <- plot_rsv_cov(
      df_avg = data_avg[[i]],
      plot_name = names(data_avg[i]),
      base_location = base_location
    )

    plots <- append(plots, list(p))
  }

  nCol <- 4
  nRow <- 2
  plots_combined <- marrangeGrob(plots, ncol = nCol, nrow = nRow)
  return(plots_combined)
}


plot_combine_flu <- function(data_avg, sample, gene, base_location) {

  plots <- list()
  for (i in seq_along(data_avg)) {
    p <- plot_flu_cov(
      df = data_avg[[i]],
      sample = sample[[i]],
      gene = gene[[i]],
      base_location = base_location
    )

    plots <- append(plots, list(p))
  }

  n <- length(plots)
  nCol <- 4
  nRow <- 2
  plots_combined <- marrangeGrob(plots, ncol = nCol, nrow = nRow)
  return(plots_combined)
}

#usage example is here:

# base_location = '~/Desktop/assemblies'
# 
# locations_aa = get_data_location(base_location, type = 'aa')
# data_aa = get_data(locations_aa)
# data_aa_avg = lapply(data_aa, FUN = average_counts, type = 'aa', by_num = 5)
# names(data_aa_avg) = str_match(locations_aa, pattern = "assemblies\\/(\\w+)\\/tables")[,2]
# plot_rsv_cov(data_aa, 'test', base_location = base_location)
# 
# file_names = tibble(
#   sample = str_match(locations_aa, pattern = "(\\w+)/tables")[, 2],
#   file_name = str_match(locations_aa, pattern = "[A|B]_\\w{2,3}.+")) %>% 
#   mutate(gene = str_replace(file_name, pattern = "-allAlleles.txt", replacement = "")) %>%
#   mutate(sample_gene = paste0(sample, "/", gene))
# 
# a = plot_combine_flu(data_avg = data_aa_avg,
#                  gene = file_names$gene,
#                  sample = file_names$sample,
#                  base_location = base_location)
# 
# ggsave(filename =  paste0("depthReport.pdf"),
#        plot = a, device = "pdf",
#        dpi = 300,  width = 420, height = 297, units = 'mm')
