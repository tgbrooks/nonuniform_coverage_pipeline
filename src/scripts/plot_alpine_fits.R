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

sample_info <- read_tsv(snakemake@input$sample_info) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

# Plot for the paper
if (snakemake@wildcards$tissue == "liver") {
    select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973")
} else if (snakemake@wildcards$tissue == "UHR") {
    select_genes <- c("ENST00000380680", "ENST00000604000", "ENST00000426077")
} else if (snakemake@wildcards$tissue == "UHR_degraded") {
    select_genes <- c("ENST00000053468", "ENST00000253063", "ENST00000272233")
} else {
    select_genes <- c()
}
if (length(select_genes) > 0) {
    print(select_genes)
    select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
        left_join(sample_info, by=join_by(sample == ID)) %>%
        group_by(sample, gene) %>%
        mutate(rel_read_depth = actual / max(actual)) %>%
        ungroup()
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
}

if (("rin_score" %in% colnames(sample_info)) & (any(~is.na(sample_info$rin_score)))) {
    # BY RIN SCORE
    select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
        left_join(sample_info, by=join_by(sample == ID)) %>%
        group_by(sample, gene) %>%
        mutate(rel_read_depth = actual / max(actual)) %>%
        ungroup()
    ggplot(
            data = select_cov,
            aes(x=pos, y=rel_read_depth, group=sample)
        ) +
        facet_grid(
            cols = vars(gene),
            scales = "free",
        ) +
        geom_path(
            aes(color = rin_score),
        ) +
        scale_color_viridis_c(option = "viridis", limits=c(0,10)) +
        labs(
            x = "Position",
            y = "Normalized Read Depth",
            color = "RIN score"
        ) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    ggsave(
        paste(outdir, "selected_genes.by_rin_score.png", sep="/"),
        width = 7,
        height = 2.5,
    )
}

# Make the plots
for (gene in unique(cov_table$gene)) {
  print(paste("Plotting", gene))
  gene_cov <- cov_table[cov_table$gene == gene,] %>%
        left_join(sample_info, by=join_by(sample == ID))
  print(gene_cov |> colnames())
  ggplot(
          data = gene_cov,
          aes(x=pos, y=actual)
      ) +
      facet_grid(
          cols=vars(sample_num),
          rows=vars(study),
          scales = "free",
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
      )
  ggsave(
      paste(outdir, paste("trace.", gene, ".png", sep=""), sep="/"),
      width = 7,
      height = 7,
  )
}
