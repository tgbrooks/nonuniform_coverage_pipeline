library(readr)
library(dplyr)
library(ggplot2)

sample_ids = c(
    "SRX4080514",
    "SRX4080520",
    #"SRX6902444",
    #"SRX6902450",
    #"SRX3304756",
    #"SRX3304764",
    "SRX4393368",
    "SRX4393369",
    "SRX16386863",
    "SRX16386864",
    "SRX14468350",
    "SRX14468347",
    "SRX13396189",
    "SRX13396186",
    "SRX11694510",
    "SRX11694499"
)
tissue <- "liver"


cov_tables_files = paste0("results/alpine_fits/", tissue, "/", sample_ids, ".coverage_table.txt")

cov_tables <- list()
for (file in cov_tables_files) {
  print(paste("Reading", file))
  cov_tables[[length(cov_tables)+1]] <- read_delim(file, delim="\t")
}
cov_table <- bind_rows(cov_tables)

sample_info <- read_tsv(paste0("results/", tissue, ".sample_info.txt")) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

# Plot for the paper
select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973")
select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth = actual / max(actual)) %>%
    arrange(study) %>%
    ungroup()
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=rel_read_depth)
    ) +
    facet_grid(
        rows=vars(study),
        cols=vars(gene),
        scales = "free",
    ) +
    coord_cartesian(xlim = c(1.5,2)) +
    geom_path(
        aes(color = sample_num),
    ) +
    labs(
        x = "Position (kb)",
        y = "Normalized read depth",
        color = "Sample\nnumber"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "zoom.liver.selected_genes.replicates.png", sep="/"),
    width = 7,
    height = 7,
)
