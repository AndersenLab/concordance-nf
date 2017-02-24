library(tidyverse)

df <- dplyr::bind_rows(lapply(dir(pattern = "*.tsv"), function(x) {
  readr::read_tsv(x, col_names = c("CHROM", "POS", "fq", "sm", "gt")) %>%
  tidyr::unite("CHROM_POS", CHROM, POS, sep = "_")
})) %>%
  dplyr::mutate(gt = as.integer(substring(gt, 1, 1)))

fq_sm <- df %>% dplyr::select(fq, sm) %>% dplyr::distinct()

d <- df %>% dplyr::select(-sm) %>%
       tidyr::spread(fq, gt) %>%
       dplyr::select(-CHROM_POS)

# Alt only
alt_only<- apply(d,2,function(x)colSums(x==d, na.rm = T))
total <- apply(d,2,function(x)colSums(!is.na(d) & !is.na(x)))

alt_only <- alt_only %>% dplyr::tbl_df() %>% 
  dplyr::mutate(b = colnames(alt_only)) %>%
  dplyr::select(b, everything()) %>%
  tidyr::gather(fq, concordant_sites, -b) %>% 
  dplyr::rename(a = fq)

total <- total %>% dplyr::tbl_df() %>% 
  dplyr::mutate(b = colnames(total)) %>%
  dplyr::select(b, everything()) %>%
  tidyr::gather(fq, total_sites, -b) %>% 
  dplyr::rename(a = fq)
  
df <- dplyr::left_join(alt_only, total) %>%
  dplyr::mutate(concordance = concordant_sites/total_sites) %>%
  dplyr::left_join(fq_sm %>% dplyr::rename(a = fq, sm1 = sm)) %>% 
  dplyr::left_join(fq_sm %>% dplyr::rename(b = fq, sm2 = sm))
