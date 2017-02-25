library(tidyverse)

gt_dict = list("0/0" = 0, "1/1" = 1)

df <- dplyr::bind_rows(lapply(dir(pattern = "rg_gt.tsv"), function(x) {
  readr::read_tsv(x, col_names = c("CHROM", "POS", "gt", "fq", "SM")) %>%
  tidyr::unite("CHROM_POS", CHROM, POS, sep = "_") %>%
  dplyr::filter(gt %in% c("0/0", "1/1")) %>%
  dplyr::mutate(gt = gt_dict[gt][[1]])
})) 

SM <- df$SM[[1]]

d <- df %>%
       tidyr::spread(fq, gt) %>%
       dplyr::select(-CHROM_POS, -SM)

if(ncol(d) == 1) {
  q()
} else {

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
  dplyr::mutate(concordance = concordant_sites/total_sites, SM = SM)

readr::write_tsv(df, "out.tsv", col_names = F)
}