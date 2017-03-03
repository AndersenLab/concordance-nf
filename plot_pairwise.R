library(tidyverse)
df <- readr::read_tsv("out.tsv", col_names = c("CHROM", "POS", "PAIR")) %>%
      dplyr::mutate(bin = round(POS/1E6)) %>%
      dplyr::group_by(CHROM,bin) %>%
      dplyr::summarize(count=n())

if (max(df$count) > 250) {
    set_max = NA
else {
    set_max = 250
}

ggplot(df, aes(x = bin, y = count)) + 
  geom_bar(aes(fill = count >= 50), stat='identity') +
  scale_fill_manual(values = c("gray", "red")) +
  facet_grid(. ~ CHROM, scales = "free_x") +
  theme_bw()  +
  scale_y_continuous(limits = c(0, set_max))

ggsave("out.png", width = 12, height = 3)