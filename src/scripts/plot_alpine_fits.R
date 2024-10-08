library(readr)
library(dplyr)
library(ggplot2)

outdir <- snakemake@output$outdir

dir.create(outdir)

cov_tables <- list()
for (file in snakemake@input$cov_tables) {
  print(paste("Reading", file))
  cov_tables[[length(cov_tables)+1]] <- read_delim(file, delim="\t")
}
cov_table <- bind_rows(cov_tables)

sample_info <- read_tsv("results/sample_info.txt") |>
    group_by(study) |>
    mutate(sample_num = row_number()) |>
    ungroup() |>
    mutate(sample_num = as.factor(sample_num))

# Plot for the paper
select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973")
select_cov <- cov_table[cov_table$gene %in% select_genes,] |>
    left_join(sample_info, by=join_by(sample == ID)) |>
    group_by(sample, gene) |>
    mutate(rel_read_depth = actual / max(actual)) |>
    ungroup()
print(select_cov)
ggplot(
        data = select_cov,
        aes(x=pos, y=rel_read_depth)
    ) +
    facet_grid(
        rows=vars(study),
        cols=vars(gene),
        scales = "free",
    ) +
    geom_path(
        aes(color = sample_num),
    ) +
    labs(
        x = "Position",
        y = "Normalized Read Depth",
        color = "Sample\nNumber"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste(outdir, "selected_genes.replicates.png", sep="/"),
    width = 7,
    height = 7,
)

# Make the plots
for (gene in unique(cov_table$gene)) {
  print(paste("Plotting", gene))
  gene_cov <- cov_table[cov_table$gene == gene,]
  ggplot(
          data = gene_cov,
          aes(x=pos, y=actual)
      ) +
      facet_wrap(
          ~sample,
          scales = "free_y",
          ncol = 1,
      ) +
      geom_path(
          color = "black",
      ) +
      geom_path(
          aes(y = predicted),
          color = "red"
      ) +
      labs(
          x = "Position",
          y = "Read Depth",
          title = "Alpine model fits"
      ) +
      theme(
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
      )
  ggsave(
      paste(outdir, paste("trace.", gene, ".png", sep=""), sep="/"),
      width = 7,
      height = 7,
  )
}
