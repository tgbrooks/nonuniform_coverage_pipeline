library(tidyverse)

outdir <- snakemake@output$outdir
dir.create(outdir)

cov_table <- read_tsv(snakemake@input$cov)

sample_info <- read_tsv(snakemake@input$sample_info) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

# Selected for the mouse liver sample set
select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973")

select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth_cov = cov / max(cov)) %>%
    mutate(rel_read_depth_cov_nonunique = cov_nonunique / max(cov_nonunique)) %>%
    mutate(rel_read_depth_cov_disconnected = cov_disconnected / max(cov_disconnected)) %>%
    ungroup()

select_cov_long <- pivot_longer(
    select_cov %>% select(
        sample,
        study,
        sample_num,,
        gene,
        pos,
        starts_with('rel_read_depth'),
    ),
    cols = starts_with('rel_read_depth'),
    names_to = "cov_type",
    names_prefix = "rel_read_depth",
    values_to = "cov",
)

ggplot(
        data = select_cov_long %>% filter(sample_num == 1),# only first samples for simplicity
        aes(x=pos, y=cov, color = cov_type)
    ) +
    facet_grid(
        cols=vars(gene),
        rows=vars(sample),
        scales = "free",
    ) +
    scale_color_viridis_d(option = "viridis", direction=-1) +
    geom_path() +
    labs(
        x = "Position",
        y = "Normalized read depth",
        color = "Coverage type",
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste(outdir, "selected_genes.by_cov_type.png", sep="/"),
    width = 8,
    height = 6,
)
