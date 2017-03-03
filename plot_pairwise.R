library(tidyverse)
df <- readr::read_tsv("out.tsv", col_names = c("CHROM", "POS", "PAIR"))

ggplot(df, aes(x = POS)) + 
  geom_histogram() +
  facet_grid(. ~ CHROM)

ggsave("out.png", width = 12, height = 3)