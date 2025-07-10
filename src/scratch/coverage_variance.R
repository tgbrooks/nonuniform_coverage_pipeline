library(readr)
library(dplyr)
library(magrittr)

cov <- read_tsv("results/transcript_coverage/liver.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/liver.sample_info.txt")

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    left_join(
        sample_info %>% select(sample = ID, study),
        by = "sample",
    )

# Position-wise fit
# These are done base-to-base since we want to account for positional variation
# but the full model is slow to run if we include a separate term for each base
# So instead we run one lm() for each base but explicitly subtract out the grand mean
# and add a dummy 'position' varaible that (since intercept is not included)
# accounts for the difference from the grand mean, thereby measuring variance of
# position.
grand_mean <- cov_downsampled$cov_normalized %>% mean()
models = list(
    base = cov_normalized - grand_mean ~ 0,
    position = cov_normalized - grand_mean ~ 1,
    study = cov_normalized - grand_mean ~ study,
    bio_rep = cov_normalized - grand_mean ~ sample # perfect fit
)
model_fit <- function(data, keys) {
    fits <- lapply(models, function(model) { lm(model, data) })
    names(fits) <- NULL
    aov <- do.call(anova, fits)
    rownames(aov) <- names(models)
    tibble(
        gene_id = keys$gene,
        pos = keys$pos,
        variable = rownames(aov),
        sum_sq = aov$`Sum of Sq`,
    )
}
squared_sums <- cov_downsampled %>% group_by(gene, pos) %>% group_map(model_fit) %>% bind_rows()

overall <- squared_sums %>%
    filter(variable != 'base') %>%
    group_by(variable) %>%
    summarize(sum_sq = sum(sum_sq)) %>%
    mutate(pct_variance = 100*sum_sq / sum(sum_sq)) %>%
    arrange(factor(variable, levels=names(models)))


write_tsv(overall, "results/scratch/coverage.GEO.percent_variance_explained.txt")


## SAME BUT FOR THE SEQC / UHR DATASET
cov <- read_tsv("results/transcript_coverage/UHR.transcript_coverage.txt.gz")
sample_info <- read_tsv("results/UHR.sample_info.txt") %>%
    mutate(
        library_prep_batch = case_when(
            library_id %in% c(1,2,3,4) ~ str_split_i(study, "_", 2),
            library_id == 5 ~ "SEQC",
        ),
        sequencing_batch = str_split_i(study, "_", 2),
    )

cov_downsampled <- cov %>%
    group_by(gene, sample) %>%
    mutate( mean_gene_cov = mean(cov) ) %>%
    filter( pos %% 100 == 50 ) %>%
    mutate( cov_normalized = cov / mean_gene_cov) %>%
    left_join(
        sample_info %>% select(sample = ID, study),
        by = "sample",
    )

# Position-wise fit
# These are done base-to-base since we want to account for positional variation
# but the full model is slow to run if we include a separate term for each base
# So instead we run one lm() for each base but explicitly subtract out the grand mean
# and add a dummy 'position' varaible that (since intercept is not included)
# accounts for the difference from the grand mean, thereby measuring variance of
# position.
grand_mean <- cov_downsampled$cov_normalized %>% mean()
models = list(
    base = cov_normalized - grand_mean ~ 0,
    position = cov_normalized - grand_mean ~ 1,
    library_prep = cov_normalized - grand_mean ~ library_prep_batch,
    sequencing_batch = cov_normalized - grand_mean ~ library_prep_batch + sequencing_batch,
    tech_rep = cov_normalized - grand_mean ~ sample
)
model_fit <- function(data, keys) {
    fits <- lapply(models, function(model) { lm(model, data) })
    names(fits) <- NULL
    aov <- do.call(anova, fits)
    rownames(aov) <- names(models)
    tibble(
        gene_id = keys$gene,
        pos = keys$pos,
        variable = rownames(aov),
        sum_sq = aov$`Sum of Sq`,
    )
}
squared_sums <- cov_downsampled %>% group_by(gene, pos) %>% group_map(model_fit) %>% bind_rows()

overall <- squared_sums %>%
    filter(variable != 'base') %>%
    group_by(variable) %>%
    summarize(sum_sq = sum(sum_sq)) %>%
    mutate(pct_variance = 100*sum_sq / sum(sum_sq)) %>%
    arrange(factor(variable, levels=names(models)))


write_tsv(overall, "results/scratch/coverage.UHR.percent_variance_explained.txt")
