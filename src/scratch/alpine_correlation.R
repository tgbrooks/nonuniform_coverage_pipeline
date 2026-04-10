library(tidyverse)

#sample_ids <- c("SRX11694499", "SRX14468347")
#tissue <- "liver"
sample_info <- read_tsv("results/liver.sample_info.txt")
sample_ids <- sample_info$ID

training_transcripts <- read_tsv('data/liver/training_set_transcripts.txt')

#training_transcripts <- read_tsv('data/liver/training_set_transcripts.txt')# LOAD DATA
temp <- list()
for (sample_id in sample_ids) {
    cov_file <- paste0("results/alpine_fits/liver/", sample_id, ".coverage_table.txt")
    if (file.exists(cov_file)) {
        temp[[length(temp)+1]] <- read_tsv(cov_file) |>
            mutate(sample_id = sample) |>
            rename(gene_id = gene, cov = actual, alpine=predicted) |>
            group_by(gene_id) |>
            mutate(mean_cov = mean(cov))
    }
}
coverage <- bind_rows(temp)


sample_ids <- sample_ids[sample_ids %in% coverage$sample_id]

### COMPUTE SMOOTED COVERAGE
HIGH_EXPRESSION_THRESHOLD <- 50


#### COMPUTE THE DISTANCE METRICS
temp <- list()
for (this_sample_id in sample_ids) {
    cov <- coverage |>
        filter(
            sample_id == this_sample_id,
            cov >= HIGH_EXPRESSION_THRESHOLD,
            !(gene_id %in% training_transcripts$transcript_id), # only use TEST transcripts
        )
    temp[[length(temp)+1]] <- tibble(
        sample_id = this_sample_id,
        cov_corr = cor.test(
            cov$cov / cov$mean_cov,
            cov$alpine / cov$mean_cov
        )$estimate,
    )
}
results <- bind_rows(temp)

print(results)

