#!/usr/bin/env Rscript

library(tidyverse)
# Read data from Rscript input
args <- commandArgs(trailingOnly=TRUE)

df_raw <- readr::read_tsv(args[3], col_names = TRUE)

df_raw_1 <- df_raw %>%
  dplyr::select(CHROM, POS, args[1], args[2]) #BRC20113_JU1530#

names(df_raw_1) <- c("CHROM", "POS", "A", "B")

df <- df_raw_1 %>%
  dplyr::mutate(bin = round(POS/1E6)) %>%
  dplyr::filter(A != "./.", B != "./.", A != "0/1", B != "0/1") %>%
  dplyr::mutate(concordant = (A == B), discordant = (A != B)) %>%
  dplyr::group_by(CHROM,bin) %>%
  dplyr::summarize(count=n(),
                   discordant_count = sum(discordant),
                   concordant = sum(concordant)/count,
                   discordant = sum(discordant)/count) %>%
  dplyr::filter(CHROM != "MtDNA")

## Plot pairwise disconcordance across Chromosomes
ggplot(df, aes(x = bin, y = discordant)) + 
    geom_bar(aes(fill = discordant >= 0.02), stat='identity') +
    scale_fill_manual(values = c("gray", "red")) +
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") +
    theme_bw() + 
    #scale_y_continuous(limits = c(0, 0.10)) +
    labs(x = "Position", y = "Discordant SNPs (%)", title = glue::glue("{args[1]} VS {args[2]}"))
  ggsave(glue::glue("{args[1]}-{args[2]}.disconcordance.png"), width = 12, height = 3)#

## Plot histogram for concordance distribution
ggplot(df) + 
    aes(x = concordant) +
    geom_histogram(aes(fill = discordant >= 0.02), binwidth = 0.01) +
    scale_fill_manual(values = c("#808080", "#0080FF")) +
    labs(x = "Condordance", y = "count", title = glue::glue("{args[1]} VS {args[2]}"))
  ggsave(glue::glue("{args[1]}-{args[2]}.hist.png"), width = 12, height = 3)
