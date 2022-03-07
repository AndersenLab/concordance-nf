# library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

df <- readr::read_tsv("out.tsv", col_names = c("CHROM", "POS", "A", "B")) %>%
  dplyr::mutate(bin = round(POS/1E6)) %>%
  dplyr::filter(A != "./.", B != "./.", A != "0/1", B != "0/1") %>%
  dplyr::mutate(concordant = (A == B), discordant = (A != B)) %>%
  dplyr::group_by(CHROM,bin) %>%
  dplyr::summarize(count=n(),
                   discordant_count = sum(discordant),
                   concordant = sum(concordant)/count,
                   discordant = sum(discordant)/count) %>%
  dplyr::filter(CHROM != "MtDNA")

ggplot(df, aes(x = bin, y = discordant)) + 
  geom_bar(aes(fill = discordant >= 0.02), stat='identity') +
  scale_fill_manual(values = c("gray", "red")) +
  facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") +
  theme_bw()  +
  scale_y_continuous(limits = c(0, 0.10)) +
  labs(x = "Position", y = "Discordant SNPs (%)") +
  ggtitle(paste0(args[2], ": ", args[1], ", match to previous isotype ", args[3]))

ggsave(paste0(args[2], ".", args[3], ".", str_replace(args[1], ",", "_"), ".png"), width = 12, height = 3)
