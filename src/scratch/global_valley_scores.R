library(tidyverse)
library(stringr)


# UHR degraded (RIN score association)
sample_info <- read_tsv("results/UHR_degraded.sample_info.txt", show_col_types=FALSE)
gvs <- read_tsv("results/UHR_degraded/corrected_global_valley_score.txt", show_col_types=FALSE) |>
    left_join(sample_info, by = join_by(sample == ID))

ggplot(gvs, aes(x=rin_score, y=global_valley_score)) +
    geom_smooth(method='lm', formula= y~x) +
    geom_point() +
    labs(
         x = "RIN score",
         y = "Global valley score",
    )
ggsave("results/scratch/UHR_degraded.global_valley_scores.png", width=3, height=3)

rin_pvalue <- summary(lm(global_valley_score ~ rin_score, data=gvs))$coefficients[2,4]
message(paste("GVS ~ RIN score pvalue", rin_pvalue))

# Testis (PCR cycle count)
gvs <- read_tsv("results/testis/global_valley_scores.txt", show_col_types=FALSE)

ggplot(gvs, aes(x=PCR_cycle_count, y=global_valley_score)) +
    geom_smooth(method='lm', formula= y~x) +
    geom_point() +
    labs(
         x = "PCR cycle number",
         y = "Global valley score",
    )
ggsave("results/scratch/testis.PCR_cycle_number.global_valley_scores.png", width=3, height=3)

PCR_cycle_pvalue <- summary(lm(global_valley_score ~ PCR_cycle_count, data=gvs))$coefficients[2,4]
message(paste("GVS ~ PCR cycle count pvalue", PCR_cycle_pvalue))


# By gene length
library(ensembldb)
liver <- read_tsv("results/liver/global_valley_scores.txt", show_col_types=FALSE)
# read in the transcript valley scores (by transcript)
temp <- list()
for (sample_id in liver$sample) {
    temp[[length(temp)+1]] <- read_tsv(paste0("data/", sample_id, "/transcript_valley_score.txt"), show_col_types=FALSE) |>
        dplyr::select(sample_id = sample, transcript_valley_score = transcript_valley_score, transcript_id = gene_id)
}
tvs <- bind_rows(temp)

ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
txdb <- EnsDb(ensembldb_sqlite)
transcript_lengths <- lengthOf(txdb, of = "tx")
tx_lens <- tibble(
    transcript_id = names(transcript_lengths),
    transcript_length = transcript_lengths,
)

tvs <- tvs |> left_join(tx_lens, by=join_by(transcript_id == transcript_id))

ggplot(tvs, aes(x=transcript_length, y=transcript_valley_score)) +
    facet_wrap(vars(sample_id)) +
    geom_point()
ggsave("results/scratch/liver.TVS_by_length.png", width=5, height=5)
txdf <- transcripts(txdb, return.type="DataFrame")
mt_genes <- (txdf |> as_tibble() |> dplyr::filter(str_detect(tx_external_name, "mt-")))$transcript_id
#tab <- table(txdf$gene_id)

# By read depth
temp <- list()
for (sample_id in liver$sample) {
    temp[[length(temp)+1]] <- read_tsv(paste0("data/", sample_id, "/transcript_coverage.txt"), show_col_types=FALSE) |>
        dplyr::select(sample_id = sample, transcript_id = gene, cov = cov)
}
cov <- bind_rows(temp)

tx_cov <- cov |> group_by(sample_id, transcript_id) |> summarize(cov = mean(cov), .groups="drop")
tvs2 <- tvs |> left_join(tx_cov, join_by(sample_id==sample_id, transcript_id == transcript_id)) |>
    mutate(
           mitochondrial = transcript_id %in% mt_genes,
           )

ggplot(tvs2, aes(x=cov, y=transcript_valley_score, color=transcript_length, shape=mitochondrial)) +
    facet_wrap(vars(sample_id), scales="free") +
    scale_color_viridis_c() +
    scale_x_log10() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_point()
ggsave("results/scratch/liver.TVS_by_read_depth.png", width=7, height=7)
