library(tidyverse)


# UHR degraded (RIN score association)
sample_info <- read_tsv("results/UHR_degraded.sample_info.txt", show_col_types=FALSE)
gvs <- read_tsv("results/UHR_degraded/corrected_global_valley_score.txt", show_col_types=FALSE) |>
    left_join(sample_info, by = join_by(sample == ID))

ggplot(gvs, aes(x=rin_score, y=global_valley_score)) +
    geom_smooth(method='lm', formula= y~x) +
    geom_point() +
    labs(
         x = "RIN score",
         y = "Global valley score",
    )
ggsave("results/scratch/UHR_degraded.global_valley_scores.png", width=3, height=3)

rin_pvalue <- summary(lm(global_valley_score ~ rin_score, data=gvs))$coefficients[2,4]
message(paste("GVS ~ RIN score pvalue", rin_pvalue))

# Testis (PCR cycle count)
gvs <- read_tsv("results/testis/global_valley_scores.txt", show_col_types=FALSE)

ggplot(gvs, aes(x=PCR_cycle_count, y=global_valley_score)) +
    geom_smooth(method='lm', formula= y~x) +
    geom_point() +
    labs(
         x = "PCR cycle number",
         y = "Global valley score",
    )
ggsave("results/scratch/testis.PCR_cycle_number.global_valley_scores.png", width=3, height=3)

PCR_cycle_pvalue <- summary(lm(global_valley_score ~ PCR_cycle_count, data=gvs))$coefficients[2,4]
message(paste("GVS ~ PCR cycle count pvalue", PCR_cycle_pvalue))
