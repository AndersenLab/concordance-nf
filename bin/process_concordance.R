#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
# library(tidyverse)

# Used in calculating isotypes
stack_list <- function(x) {
    if (is.null(names(x)) == T) {
        names(x) <- as.character(1:length(x))
    }
    stack(x)
}

args <- commandArgs(trailingOnly=TRUE)

# load bam coverage data
SM_coverage <- read.delim(args[1], stringsAsFactors=FALSE) # %>% rename(strain=Sample)

# SM_coverage$coverage <- rowMeans(SM_coverage[,c("I", "II", "III", "IV", "V")], na.rm=TRUE)

SM_coverage <- select(SM_coverage, strain, coverage)

# load existing WI isotype assignment
WI <- readr::read_tsv(args[2]) 

WI %>% readr::write_tsv("WI_metadata.tsv")

existing_WI <- WI %>%
    dplyr::select(strain, isotype, latitude, longitude) %>%
    dplyr::filter(isotype != "NA")

# check if any strains have an isotype existing (for new species)
if(nrow(existing_WI) == 0) {
    existing_WI <- WI %>%
        dplyr::select(strain, isotype, latitude, longitude)
}

# define concordance cutoff for strains that belong to the same isotype
cutoff <- as.numeric(args[3])

# compute pairwise concordance and whether isotype
gtcheck <- readr::read_tsv("gtcheck.tsv") %>%
    dplyr::mutate(concordance = 1-(discordance/sites)) %>%
    dplyr::mutate(isotype = concordance > cutoff) 

# Generate strains that do not group with any other strains (single strains)
single_strains <- gtcheck %>%
    dplyr::select(i, j, isotype) %>%
    tidyr::gather(col, strain, -isotype) %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(single_strain = sum(isotype, na.rm = T)) %>% 
    dplyr::filter(single_strain==0)


single_strains <- as.data.frame(list(i = single_strains$strain,
                                     j = single_strains$strain,
                                     discordance = 0,
                                     concordance = 1,
                                     isotype = TRUE))

gtcheck <- dplyr::bind_rows(gtcheck, single_strains)

# Filter for strain pairs that are isotype
iso_gtcheck <- gtcheck %>% dplyr::filter(isotype == T)

# Generate complete strain list
strain_list <- sort(unique(c(gtcheck$i, gtcheck$j)))

# Generate isotype groups:
# For each strain, put all other strains isotype with it into its group
isotype_groups <- stack_list(unique(lapply(strain_list, function(x) {
    grouped_strains <- dplyr::filter(iso_gtcheck, (i == x | j == x)) %>%
        dplyr::select(i, j)
    sort(unique(unlist(grouped_strains)))
}))) %>%
    dplyr::rename(strain = values, group = ind) %>%
    dplyr::select(group, everything())

# Add info about coverage and collection locations
# Add column "strain_in_multiple_isotypes" to indicate if 1 strain fall in multiple groups
# Add column "unique_isotypes_per_group" to count for each new group, how many exisiting isotypes it contains
# Add column "unique_groups_per_isotype" to count for each existing isotype how many new groups it contains (note those with no existing isotypes will have NA and be counted as the same isotype)
isotype_groups <- dplyr::left_join(isotype_groups, existing_WI, by = c("strain")) %>%
    dplyr::left_join(SM_coverage) %>%
    dplyr::mutate(group = as.integer(group)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(avg_lat = abs(as.numeric(latitude) - mean(as.numeric(latitude), na.rm = T)), 
                  avg_lon = abs(as.numeric(longitude) - mean(as.numeric(longitude), na.rm = T))) %>%
    dplyr::mutate(unique_isotypes_per_group = length(unique(purrr::discard(isotype, is.na)))) %>%
    dplyr::group_by(isotype) %>%
    dplyr::mutate(unique_groups_per_isotype = length(unique(group))) %>%
    dplyr::group_by(strain) %>%
    dplyr::mutate(strain_in_multiple_isotypes = length(strain) > 1) %>%
    dplyr::mutate(location_issue = (avg_lat > 5 | avg_lon > 5)) %>%
    dplyr::select(-avg_lat,-avg_lon) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(strain_conflict = any(unique_isotypes_per_group > 1,
                                        unique_groups_per_isotype > 1,
                                        location_issue,
                                        na.rm = T))


# Save text files
readr::write_tsv(isotype_groups, paste0("isotype_groups.tsv"))

# for pairwise - only isotype groups with new strains
groups <- isotype_groups %>%
    dplyr::filter(is.na(isotype))

new <- isotype_groups %>%
    dplyr::filter(group %in% unique(groups$group))
    
readr::write_tsv(new, paste0("isotype_groups_new.tsv"))


# Output problem strains
pr_strains <- unique((isotype_groups %>% 
                          dplyr::filter(strain_conflict))$strain)

pr_strain_comparison <- gtcheck %>% dplyr::filter(i %in% pr_strains, j %in% pr_strains) %>% 
    dplyr::filter(discordance < 10000) %>% 
    dplyr::left_join(existing_WI %>% dplyr::rename(iso_i = isotype) %>% dplyr::select(-latitude, -longitude), by = c("i"= "strain")) %>%
    dplyr::left_join(existing_WI %>% dplyr::rename(iso_j = isotype) %>% dplyr::select(-latitude, -longitude), by = c("j"= "strain")) %>%
    dplyr::filter(!is.na(sites))

readr::write_tsv(pr_strain_comparison, paste0("problem_strains.tsv"))


# Count total number of new isotype groups
isotypes <- length(unique(isotype_groups$group))

f <- file("isotype_count.txt")
writeLines(paste0(isotypes, sep = "\n"), f)


# Make plot of distribution of concordance
ggplot(gtcheck) +
    geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.00025) +
    scale_fill_manual(values = c("#808080", "#0080FF"))

ggsave("concordance.pdf", width = 5, height = 5)
ggsave("concordance.png", width = 5, height = 5)

ggplot(gtcheck) +
    geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.000025) +
    scale_fill_manual(values = c("#808080", "#0080FF")) +
    scale_x_continuous(limits = c(0.99, 1.0)) +
    labs(x = "Concordance", y = "Number of Comparisons") +
    geom_vline(aes(xintercept = cutoff), color = "red") +
    theme(axis.title = ggplot2::element_text(size=14, face="bold", color="black", vjust=5))

ggsave("xconcordance.pdf", width = 5, height = 5)
ggsave("xconcordance.png", width = 5, height = 5)

