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
