library(tidyverse)

cor <- read_tsv("results/UHR/coverage_distance.txt")
sample_info <- read_tsv("results/UHR.sample_info.txt")
data <- cor |>
    left_join(sample_info |> select(sample_id1 = ID, study1 = study), by = "sample_id1") |>
    left_join(sample_info |> select(sample_id2 = ID, study2 = study), by = "sample_id2") |>
    mutate(
        study1 = str_extract(study1, "PRJNA208369_([A-Z]+)", group=1),
        study2 = str_extract(study2, "PRJNA208369_([A-Z]+)", group=1),
    )

within_corr <- data |> filter(study1 == study2) |> group_by(study1) |> summarize(cov_corr = mean(cov_corr), local_corr = mean(local_corr))
out_corr <- data |> filter(study1 != study2) |> group_by(study1, study2) |> summarize(cov_corr = mean(cov_corr), local_corr = mean(local_corr), .groups="drop")

ggplot(data, aes(x = study2, y = local_corr, color = study2)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        data = data |> filter(study1 != study2),
        position=position_jitter(width=0.3, height=0),
        shape = "circle open"
    ) +
    geom_jitter(
        data = data |> filter(study1 == study2),
        position=position_jitter(width=0.3, height=0),
    ) +
    labs(
         color = "Site",
         x = "Site",
         y = "Local correlation"
    ) +
    ylim(0.7, 1.0)
ggsave("results/scratch/UHR.coverage_distance.png", width=5, height=5)


cor <- read_tsv("results/liver/coverage_distance.txt")
sample_info <- read_tsv("results/liver.sample_info.txt")
data <- cor |>
    left_join(sample_info |> select(sample_id1 = ID, study1 = study), by = "sample_id1") |>
    left_join(sample_info |> select(sample_id2 = ID, study2 = study), by = "sample_id2")

within_corr <- data |> filter(study1 == study2) |> group_by(study1) |> summarize(cov_corr = mean(cov_corr), local_corr = mean(local_corr))
out_corr <- data |> filter(study1 != study2) |> group_by(study1, study2) |> summarize(cov_corr = mean(cov_corr), local_corr = mean(local_corr), .groups="drop")


ggplot(data, aes(x = study2, y = cov_corr, color = study2)) +
    facet_grid(
        rows = "study1",
    ) +
    geom_jitter(
        data = data |> filter(study1 != study2),
        position=position_jitter(width=0.2, height=0),
        shape = "circle open"
    ) +
    geom_jitter(
        data = data |> filter(study1 == study2),
        position=position_jitter(width=0.2, height=0),
    ) +
    labs(
         color = "Site",
         x = "Site",
         y = "Global correlation"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(0.6, 1.0)
ggsave("results/scratch/liver.coverage_distance.png", width=5, height=7)
