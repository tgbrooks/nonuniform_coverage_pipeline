### SCRIPT TO ASSESS WHETHER EXONS HAVE AN EFFECT ON COVERAGE ###
### here we model expression as a linear model with exon ID (as a factor)
### so that each exon can have a different value. This is compared to models with
### randomly permuted exon boundaries.
library(ensembldb)
library(tidyverse)


sample_ids <- snakemake@params$sample_ids
ensembldb_sqlite <- snakemake@input$ensdb
outdir <- snakemake@output$outdir
cov_files <- snakemake@input$cov_tables
transcript_file <- snakemake@input$transcripts
sample_info_file <- snakemake@input$sample_info
tissue <- snakemake@wildcards$tissue

set.seed(1)

## TESTING VALUES
#sample_ids <- c("SRX4080514", "SRX4080520", "SRX13396186", "SRX13396189") 
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#transcript_file <- "data/liver/high_expressed_single_isoform_genes.txt"
#outdir <- "results/temp/"
#cov_files <- paste0("results/alpine_fits/liver/", sample_ids, ".coverage_table.txt")
#sample_info_file <- "results/liver.sample_info.txt"
#tissue <- "liver"
#outdir <- "results/temp"
##########

dir.create(outdir)

txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

ebt0 <- exonsBy(txdb, by="tx")


dir.create(outdir)

cov_tables <- list()
for (file in cov_files) {
  print(paste("Reading", file))
  cov_tables[[length(cov_tables)+1]] <- read_delim(file, delim="\t")
}
cov_table <- bind_rows(cov_tables)

high_exp_genes <- read_tsv(transcript_file)

sample_info <- read_tsv(sample_info_file) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

genome_to_transcript <- function(genome_pos, exons, strand) {
    out <- rep(NA, length(genome_pos))
    for (i in seq_along(genome_pos)) {
        pos <- genome_pos[[i]]
        offsets <- exons |> arrange(exon_start) |> mutate(offset = c(0, cumsum(length)) |> head(-1))
        overlap <- offsets |> filter(exon_start <= pos, exon_end >= pos)
        value <- if (nrow(overlap) > 0) {
            pos - overlap$exon_start[[1]] + overlap$offset[[1]] 
        } else {
            NA
        }
        if (strand == 2) {
            # Reverse strand, we start from the other end
            value <- sum(exons$length) - value
        }
        out[[i]] <- value
    }
    return(out)
}

get_exon_cds_info <- function(tx_id) {
    this_tx_id <- tx_id
    tx <- ebt0[[tx_id]]
    tx_contig <- (seqnames(tx) |> as.character())[1]
    tx_strand <- (strand(tx) |> as.vector() |> case_match( "+" ~ 1, "-"  ~ 2))[1]
    exon_order <- if (tx_strand == 1) {1} else {-1}
    exons <- tibble(
        exon_start = start(tx),
        exon_end = end(tx),
        exon_id = tx$exon_id,
        length = lengths(tx),
    )
    this_txdf <- txdf |> as_tibble() |> filter(tx_id == this_tx_id)

    # Take CDS and exon boundaries in transcript coordinates
    # accounting for strand
    cds_start1 <- this_txdf$tx_cds_seq_start[[1]]
    cds_end1 <- this_txdf$tx_cds_seq_end[[1]]
    cds_start2 <- genome_to_transcript(cds_start1, exons, tx_strand)
    cds_end2 <- genome_to_transcript(cds_end1, exons, tx_strand)
    cds_start <- min(cds_start2, cds_end2)
    cds_end <- max(cds_start2, cds_end2)
    exon_start1 = genome_to_transcript(exons$exon_start, exons, tx_strand)
    exon_end1 = genome_to_transcript(exons$exon_end, exons, tx_strand)
    tx_exons <- exons |> mutate(
        exon_start = pmin(exon_start1, exon_end1),
        exon_end = pmax(exon_start1, exon_end1),
    ) |>
        mutate(parity = row_number() %% 2)

    # Split exons where the CDS starts
    exons2 <- bind_rows(
            tx_exons |> filter(exon_start < cds_start) |> mutate(exon_end = min(exon_end, cds_start), type="utr"),
            tx_exons |> filter(exon_end >= cds_start, exon_start <= cds_end) |> mutate(exon_start = pmax(exon_start, cds_start), exon_end = pmin(exon_end, cds_end), type="cds"),
            tx_exons |> filter(exon_end > cds_end) |> mutate(exon_start = pmax(exon_start, cds_end), type="utr"),
        )
    return(exons2)
}

# Get all exon boundaries
temp <- list()
for (tx_id in high_exp_genes$transcript_id) {
    temp[[length(temp)+1]] <- get_exon_cds_info(tx_id) |> mutate(tx_id = tx_id)
}
exon_info <- bind_rows(temp) |>
    group_by(tx_id) |> 
    mutate(exon_number = row_number()) |>
    ungroup()

compute_exon_lm <- function(exon_info) {
    # Compute linear model `cov ~ exon` for each transcript and sample
    # returns the explained sum-squares of the regression

    # Expand exons out into base positions
    # so this maps base position (in transcriptome) to the exon number
    bp_exon_info <- exon_info |>
        group_by(tx_id, exon_number, type, exon_id, exon_start, exon_end) |>
        reframe(
            pos = seq(exon_start, exon_end)
        )

    # Select subset of positions for modelling at for efficiency
    chosen_positions <- bp_exon_info |>
        group_by(tx_id, exon_number) |>
        sample_n(10, replace=TRUE) |>
        distinct() |>
        ungroup()

    # Join positions with coverage
    model_data <- chosen_positions |>
        select(tx_id, exon_id, pos) |>
        left_join(cov_table, join_by(tx_id == gene, pos == pos), relationship="many-to-many") |>
        drop_na(sample)

    # Perform linear models
    model_expr <- function(data) {
        if (length(data$exon_id |> unique()) == 1) {
            return(NA)
        }
        res <- lm(actual ~ as.factor(exon_id), data)
        return(anova(res)[['Sum Sq']][[1]])
    }
    lm_res <- model_data |>
        group_by(sample, tx_id) |>
        summarize(
            sum_sq = model_expr(cur_data()),
            .groups = "drop"
        )
    return(lm_res)
}

# Run on the original data (true exon boundaries)
lm_res <- compute_exon_lm(exon_info)

# redo this with permuted exons
temp <- list()
for (i in seq(100)) {
    perm_exon_info <- exon_info |>
        mutate(len = exon_end - exon_start + 1) |>
        group_by(tx_id) |>
        mutate(
            new_len = sample(len, length(len)),
        ) |>
        mutate(
            exon_start = c(1, cumsum(new_len)) |> head(-1),
            exon_end = cumsum(new_len)
        )
    lm_res_ <- compute_exon_lm(perm_exon_info)
    temp[[length(temp) + 1]] <- lm_res_ |>
        mutate(
               perm_num = i,
            )
}
perm_lm_res <- bind_rows(temp)


# Generate a comparison figure
data <- perm_lm_res |>
    left_join(
        lm_res |>
            select(
                tx_id = tx_id,
                sample = sample,
                orig_sum_sq = sum_sq,
            ),
        join_by(tx_id == tx_id, sample == sample)
    ) |>
    mutate(
        norm_sum_sq = sum_sq / orig_sum_sq
    )
ggplot(data) +
    facet_wrap(vars(sample)) +
    geom_point(
        aes(x=tx_id, y=norm_sum_sq),
        alpha=0.1
    ) +
    geom_hline(
        yintercept = 1,
        color = "red"
    ) +
    scale_y_log10() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0(outdir, '/exon_effect.png'))

# Simpler summary
# Take the permutation p-value of each gene-sample and plot those
perm_p_values <- data |>
    group_by(sample, tx_id) |>
    summarize(
        perm_p = (sum(sum_sq > orig_sum_sq) + 1) / (length(sum_sq)+1),
        .groups="drop"
    )
ggplot(perm_p_values) +
    facet_wrap(vars(sample)) +
    geom_histogram(
        aes(perm_p),
        breaks=seq(0,1,by=0.1)
    ) +
    labs(x = "permutation p-value") +
    ggtitle("Permutation p-values per gene and sample")
ggsave(paste0(outdir, '/exon_effect_perm_ps.png'))
