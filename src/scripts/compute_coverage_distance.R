library(tidyverse)

#sample_ids <- c("SRX11694499", "SRX14468347")
#tissue <- "liver"
tissue <- snakemake@wildcards$tissue
sample_ids <- snakemake@params$sample_ids

## LOAD DATA
temp <- list()
for (sample_id in sample_ids) {
    cov_file <- paste0("results/alpine_fits/", tissue, "/", sample_id, ".coverage_table.txt")
    temp[[length(temp)+1]] <- read_tsv(cov_file) |>
        mutate(sample_id = sample) |>
        rename(gene_id = gene, cov = actual)
}
coverage <- bind_rows(temp)
print(coverage)


### COMPUTE SMOOTED COVERAGE
KERNEL_SIZE <- 101
CONV_SD <- 30
HIGH_EXPRESSION_THRESHOLD <- 200

create_gaussian_kernel <- function(conv_sd, kernel_size) {
  x <- seq(-kernel_size / 2, kernel_size / 2, length.out = kernel_size)
  kernel <- exp(-0.5 * (x / conv_sd)^2)
  kernel <- kernel / sum(kernel)  # Normalize the kernel
  return(kernel)
}

get_conv_result <- function(x, conv_sd=sd) {
  conv_result <- convolve(x, create_gaussian_kernel(conv_sd, KERNEL_SIZE), type = "open")
  result_length <- length(conv_result)
  rem <- floor(KERNEL_SIZE/2)
  output <- conv_result[(rem + 1) : (result_length - rem)]
  return(output)
}

coverage <- coverage |>
    group_by(sample_id, gene_id) |>
    mutate(smoothed = get_conv_result(cov,  CONV_SD)) |>
    ungroup()


#### COMPUTE THE DISTANCE METRICS
temp <- list()
for (sample_id1 in sample_ids) {
    cov1 <- coverage |> filter(sample_id == sample_id1)
    for (sample_id2 in sample_ids) {
        if (sample_id1 == sample_id2) { next }
        cov2 <- coverage |> filter(sample_id == sample_id2)
        both <- left_join(cov1, cov2, c("gene_id", "pos"))
        both_big <- both |> filter(
            cov.x >= HIGH_EXPRESSION_THRESHOLD,
            cov.y >= HIGH_EXPRESSION_THRESHOLD
        )
        temp[[length(temp)+1]] <- tibble(
            sample_id1 = sample_id1,
            sample_id2 = sample_id2,
            cov_corr = cor.test(both_big$cov.x, both_big$cov.y)$estimate,
            local_corr = cor.test(both_big$cov.x - both_big$smoothed.x, both_big$cov.y - both_big$smoothed.y)$estimate,
        )
    }
}
results <- bind_rows(temp)

results |> write_tsv(paste0("results/", tissue, "/coverage_distance.txt"))
