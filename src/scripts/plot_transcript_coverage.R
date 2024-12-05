library(readr)
library(dplyr)
library(ggplot2)

outdir <- snakemake@output$outdir
dir.create(outdir)

cov_table <- read_tsv(snakemake@input$cov)

sample_info <- read_tsv(snakemake@input$sample_info) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

select_genes <- c("ENST00000246006", "ENST00000339697", "ENST00000282388")

select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth = cov / max(cov)) %>%
    ungroup()

ggplot(
        data = select_cov,
        aes(x=pos, y=rel_read_depth)
    ) +
    facet_grid(
        cols=vars(gene),
        rows=vars(individual),
        scales = "free",
    ) +
    scale_color_viridis_c(option = "viridis", direction=-1) +
    geom_path(
        aes(color = hours),
    ) +
    labs(
        x = "Position",
        y = "Normalized read depth",
        color = "Degradation\ntime (h)",
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste(outdir, "selected_genes.by_degradation_time.png", sep="/"),
    width = 8,
    height = 5,
)
