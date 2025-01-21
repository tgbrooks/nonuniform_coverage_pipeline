library(tidyverse)

studies <- c(
    "GSE77221",
    "GSE117134",
    "GSE208768",
    "PRJNA816471",
    "PRJNA788430",
    "PRJNA753198"
)

cov_table <- read_tsv("results/transcript_coverage/liver.transcript_coverage.txt.gz")

high_exp_genes <- read_tsv("data/liver/high_expressed_single_isoform_genes.txt")

sample_info <- read_tsv("results/liver.sample_info.txt") %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))
sample_ids <- (sample_info |> filter(study %in% studies))$ID

select_genes <- c(
    "ENSMUST00000068593",
    "ENSMUST00000081777",
    "ENSMUST00000074051"
    #"ENSMUST00000082402",
    #"ENSMUST00000078869",
)

select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    filter(study %in% studies) %>%
    mutate(
        cov_normalized = cov / max(cov_nonunique),
        cov_nu_normalized = cov_nonunique / max(cov_nonunique),
    ) %>%
    ungroup() %>%
    pivot_longer(
        c(cov_normalized, cov_nu_normalized),
        names_to = "cov_type",
        values_to = "cov_normalized",
    ) |>
    mutate(
        cov_type = case_when(
            cov_type == "cov_normalized" ~ "unique only",
            cov_type == "cov_nu_normalized" ~ "any",
        )
    )

ggplot(
        data = select_cov |> filter(sample_num == 1),
        aes(x=pos, y=cov_normalized)
    ) +
    facet_grid(
        cols=vars(gene),
        rows=vars(study),
        scales = "free",
    ) +
    scale_color_viridis_c(option = "viridis", direction=-1) +
    geom_path(
              aes(color = cov_type)
    ) +
    labs(
        x = "Position",
        y = "Normalized read depth",
        color = "Coverage type",
    ) +
    scale_color_manual(values=c("black", "red")) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "multimappers.png", sep="/"),
    width = 7,
    height = 7,
)

select_genes <- c(
    "ENSMUST00000068593",
    "ENSMUST00000081777",
    "ENSMUST00000074051"
    #"ENSMUST00000082402",
    #"ENSMUST00000078869",
)

select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == ID)) %>%
    group_by(sample, gene) %>%
    filter(study %in% studies) %>%
    mutate(
        cov_normalized = cov / max(cov_nonunique),
        cov_nu_normalized = cov_nonunique / max(cov_nonunique),
    ) %>%
    ungroup() %>%
    pivot_longer(
        c(cov_normalized, cov_nu_normalized),
        names_to = "cov_type",
        values_to = "cov_normalized",
    ) |>
    mutate(
        cov_type = case_when(
            cov_type == "cov_normalized" ~ "unique only",
            cov_type == "cov_nu_normalized" ~ "any",
        )
    ) %>%
    left_join(high_exp_genes, by=join_by(gene == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", gene))

ggplot(
        data = select_cov |> filter(sample_num == 1),
        aes(x=pos, y=cov_normalized)
    ) +
    facet_grid(
        cols=vars(full_gene),
        rows=vars(study),
        scales = "free",
    ) +
    scale_color_viridis_c(option = "viridis", direction=-1) +
    geom_path(
              aes(color = cov_type)
    ) +
    labs(
        x = "Position",
        y = "Normalized read depth",
        color = "Coverage type",
    ) +
    scale_color_manual(values=c("black", "red")) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "multimappers.png", sep="/"),
    width = 7,
    height = 7,
)


# Print out some stats
ct <- cov_table |> filter(sample %in% sample_ids)
by_gene = ct |> group_by(gene) |> summarize(max((cov_nonunique+100) / (cov+100)))
