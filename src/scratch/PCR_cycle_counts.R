library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(tibble)


## Study that varied PCR cycle counts in UHR
cov <- read_tsv("results/transcript_coverage/UHR_PCR.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/UHR_PCR.sample_info.txt") |>
    mutate(
        final_amount = input_material_ng * 2^PCR_cycle_count
    )

# We don't use all genes because most have terrible coverage
# these are JUST mitochondrial genes, the only ones with high expression in every sample
selected_genes <- c(
    "ENST00000361789",
    "ENST00000361739",
    "ENST00000361381",
    "ENST00000331825",
    "ENST00000362079",
    "ENST00000361899",
    "ENST00000361624"
)

cov_downsampled <- cov %>%
    filter(gene %in% selected_genes, pos %% 10 == 5) %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov), gene_length=max(pos) ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    mutate( pos_gene = paste0(gene, pos)) %>%
    left_join(
        sample_info %>% select(sample = ID, study, PCR_cycle_count, input_material_ng, final_amount),
        by = "sample",
    )

res0 <- lm(cov_normalized ~ pos_gene, cov_downsampled)
res1 <- lm(cov_normalized ~ PCR_cycle_count * pos_gene, cov_downsampled)
res2 <- lm(cov_normalized ~ input_material_ng * pos_gene, cov_downsampled)
res3 <- lm(cov_normalized ~ (PCR_cycle_count + input_material_ng) * pos_gene, cov_downsampled)
res4 <- lm(cov_normalized ~ final_amount * pos_gene, cov_downsampled)

BIC(res0, res1, res2, res3, res4)
