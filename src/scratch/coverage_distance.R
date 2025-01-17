library(tidyverse)

shapes <- function(x) {
    c("circle", "circle open")
}

cor <- read_tsv("results/UHR/coverage_distance.txt")
sample_info <- read_tsv("results/UHR.sample_info.txt")
data <- cor |>
    left_join(sample_info |> select(sample_id1 = ID, study1 = study), by = "sample_id1") |>
    left_join(sample_info |> select(sample_id2 = ID, study2 = study), by = "sample_id2") |>
    mutate(
        study1 = str_extract(study1, "PRJNA208369_([A-Z]+)", group=1),
        study2 = str_extract(study2, "PRJNA208369_([A-Z]+)", group=1),
    ) |>
    filter(!((study1 == study2) & (sample_id1 > sample_id2))) |>
    mutate(within = case_when(study1 == study2 ~ "within", study1 != study2 ~ "between"))

### LOCAL CORRELATION UHR
ggplot(data, aes(x = study2, y = local_corr, color = study2, shape = within)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        position=position_jitter(width=0.3, height=0),
    ) +
    scale_shape_manual(values = c(1,16), name = "Correlation sites") +
    labs(
         color = "Site",
         x = "Site",
         y = "Local coverage correlation"
    ) +
    guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
    #ylim(0.7, 1.0)
ggsave("results/scratch/UHR.local.coverage_distance.png", width=5, height=5)

### GLOBAL CORRELATION UHR
ggplot(data, aes(x = study2, y = cov_corr, color = study2, shape = within)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        position=position_jitter(width=0.3, height=0),
    ) +
    scale_shape_manual(values = c(1,16), name = "Correlation sites") +
    labs(
         color = "Site",
         x = "Site",
         y = "Coverage correlation"
    ) +
    guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))
    #ylim(0.7, 1.0)
ggsave("results/scratch/UHR.coverage_distance.png", width=5, height=5)


cor <- read_tsv("results/liver/coverage_distance.txt")
sample_info <- read_tsv("results/liver.sample_info.txt")
data <- cor |>
    left_join(sample_info |> select(sample_id1 = ID, study1 = study), by = "sample_id1") |>
    left_join(sample_info |> select(sample_id2 = ID, study2 = study), by = "sample_id2") |>
    filter(!((study1 == study2) & (sample_id1 > sample_id2))) |>
    mutate(within = case_when(study1 == study2 ~ "within", study1 != study2 ~ "between"))

#### LOCAL CORRELATION LIVER
ggplot(data, aes(x = study2, y = local_corr, color = study2, shape = within)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        position=position_jitter(width=0.2, height=0),
    ) +
    scale_shape_manual(values = c(1,16), name = "Correlation studies") +
    labs(
         color = "Study",
         x = "Study",
         y = "Local coverage correlation"
    ) +
    guides(color = guide_legend(order = 2), shape = guide_legend(order = 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    #ylim(0.6, 1.0)
ggsave("results/scratch/liver.local.coverage_distance.png", width=5, height=7)

#### COVERAGE CORRELATION LIVER
ggplot(data, aes(x = study2, y = cov_corr, color = study2, shape = within)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        position=position_jitter(width=0.2, height=0),
    ) +
    scale_shape_manual(values = c(1,16), name = "Correlation studies") +
    labs(
         color = "Study",
         x = "Study",
         y = "Coverage correlation"
    ) +
    guides(color = guide_legend(order = 2), shape = guide_legend(order = 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    #ylim(0.6, 1.0)
ggsave("results/scratch/liver.coverage_distance.png", width=5, height=7)
