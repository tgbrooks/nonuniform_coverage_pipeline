library(readr)
library(dplyr)
library(ggplot2)

sample_ids = c(
    "SRX302130",
    "SRX302146",
    "SRX302162",
    "SRX302178",
    "SRX302514",
    "SRX302529",
    "SRX302544",
    "SRX302559",
    "SRX302874",
    "SRX302890",
    "SRX302906",
    "SRX302922",
    "SRX302574",
    "SRX302194",
    "SRX302938"
)
tissue <- "UHR"

high_exp_genes <- read_tsv("data/UHR/high_expressed_single_isoform_genes.txt")

### BEGIN Functions
## Convolution functions
KERNEL_SIZE <- 101
CONV_SD <- 30
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
####


cov_tables_files = paste0("results/alpine_fits/", tissue, "/", sample_ids, ".coverage_table.txt")

cov_tables <- list()
for (file in cov_tables_files) {
  print(paste("Reading", file))
  cov_tables[[length(cov_tables)+1]] <- read_delim(file, delim="\t")
}
cov_table <- bind_rows(cov_tables)

cov_table <- cov_table |>
    group_by(gene, sample) |>
    mutate(smoothed = get_conv_result(actual, conv_sd = CONV_SD))

sample_info <- read_tsv(paste0("results/", tissue, ".sample_info.txt")) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

# Plot for the paper
select_genes <- c("ENST00000380680", "ENST00000604000", "ENST00000426077")
select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth = (actual - smoothed) / max(actual - smoothed)) %>%
    left_join(high_exp_genes, by=join_by(gene == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", gene)) %>%
    arrange(study) %>%
    ungroup()
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=rel_read_depth)
    ) +
    facet_grid(
        rows=vars(study),
        cols=vars(full_gene),
        scales = "free",
    ) +
    coord_cartesian(xlim = c(1.7,2)) +
    geom_path(
        aes(color = sample_num),
    ) +
    labs(
        x = "Position (kb)",
        y = "Normalized small variations",
        color = "Sample\nnumber"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "zoom.UHR.selected_genes.replicates.png", sep="/"),
    width = 7,
    height = 7,
)
