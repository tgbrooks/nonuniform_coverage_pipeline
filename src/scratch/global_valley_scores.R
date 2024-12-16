library(tidyverse)


gvs <- read_tsv("results/UHR_degraded/global_valley_scores.txt")

ggplot(gvs, aes(x=rin_score, y=global_valley_score)) +
    geom_point() +
    labs(
         x = "RIN score",
         y = "Global valley score",
    )
ggsave("results/scratch/UHR_degraded.global_valley_scores.png", width=3, height=3)
