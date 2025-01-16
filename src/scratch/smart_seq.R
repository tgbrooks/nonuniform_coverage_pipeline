library(tidyverse)

outdir <- 'results/scratch/'

cov_table <- read_tsv('results/transcript_coverage/smart_seq.transcript_coverage.txt.gz')

sample_info <- read_tsv('results/smart_seq.sample_info.txt') %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

select_genes <- c("ENSMUST00000102476", "ENSMUST00000026565", "ENSMUST00000082421")

select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth = cov / max(cov)) %>%
    ungroup()

ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=rel_read_depth, color=seq_type, group=sample)
    ) +
    facet_grid(
        cols=vars(gene),
        scales = "free",
    ) +
    geom_path() +
    labs(
        x = "Position (kb)",
        y = "Normalized read depth",
        color = "Type",
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste(outdir, "selected_genes.smart_seq.png", sep="/"),
    width = 7,
    height = 2.5,
)
