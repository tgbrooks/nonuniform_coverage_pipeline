library(BSgenome)
library(ensembldb)
library(stringr)
library(dplyr)
library(readr)
library(tibble)

sample_ids <- snakemake@params$sample_ids
transcript_file <- snakemake@input$transcripts
ensembldb_sqlite <- snakemake@input$ensdb
coverage_file <- snakemake@input$coverage
output <- snakemake@output$outfile

# TESTING VALUES ########
#sample_ids <- c( "SRX4080514", "SRX4080520", "SRX4393368", "SRX4393369", "SRX16386863", "SRX16386864", "SRX14468350", "SRX14468347", "SRX13396189", "SRX13396186", "SRX11694510", "SRX11694499")
#transcript_file <- "data/liver/high_expressed_single_isoform_genes.txt"
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#coverage_file <- "results/transcript_coverage/liver.transcript_coverage.txt.gz"
#output <- "temp.unannotated_splicing.txt"
##########################

junction_files <- lapply(sample_ids, function(s) {paste0("data/", s, "/bam/SJ.out.tab")})

txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

selected.genes.and.transcripts <- read.csv(transcript_file, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]

library(GenomicRanges)

# For each transcript compute the mean coverage from its coverage file
coverage <- read_tsv(coverage_file)
mean_cov <- coverage |> group_by(gene, sample) |> summarize(mean_cov = mean(cov), .groups="drop")

# Load the junction data as determined by STAR
temp <- list()
for (i in seq(length(sample_ids))) {
    sample_id <- sample_ids[i]
    junction_file <- junction_files[[i]]
    col_names <- c("contig", "first_base", "last_base", "strand", "intron_motif", "annotated", "num_unique_reads_spanning", "num_multimapper_reads_spanning", "max_overhang")
    temp[[length(temp) + 1]] <- read_tsv(junction_file, col_names=col_names) |> mutate(sample_id = sample_id)
}
all_junctions <- bind_rows(temp)

# Summarize unannotated junctions over all the samples
novel_junctions <- all_junctions |>
    filter(annotated == 0) |>
    group_by(contig, first_base, last_base, strand, annotated, intron_motif, sample_id) |>
    summarize(
        max_num_unique_reads_spanning = max(num_unique_reads_spanning),
        max_num_multimapper_reads_spanning = max(num_multimapper_reads_spanning),
        max_max_overhang = max(max_overhang),
        .groups="drop",
    )

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

# For each selected transcript, check if there are any unannotated junctions within it
temp <- list()
for (tx_id in selected.txs) {
    tx <- ebt.fit[[tx_id]]
    tx_contig <- (seqnames(tx) |> as.character())[1]
    tx_strand <- (strand(tx) |> as.vector() |> case_match( "+" ~ 1, "-"  ~ 2))[1]
    exon_order <- if (tx_strand == 1) {1} else {-1}
    exons <- tibble(
        exon_start = start(tx),
        exon_end = end(tx),
        length = lengths(tx),
        exon_id = tx$exon_id,
    )
    j_start_in <- novel_junctions |>
        filter(contig == tx_contig, strand == tx_strand) |>
        inner_join(
            exons,
            by = join_by(
                    first_base >= exon_start,
                    first_base <= exon_end,
                )
            )
    j_end_in <- novel_junctions |>
        filter(contig == tx_contig, strand == tx_strand) |>
        inner_join(
            exons,
            by = join_by(
                    last_base >= exon_start,
                    last_base <= exon_end,
                )
            )
    tx_novel_j <- bind_rows(j_start_in, j_end_in) |>
        mutate(
            tx_id = tx_id,
            first_base_tx = genome_to_transcript(first_base, exons,tx_strand),
            last_base_tx = genome_to_transcript(last_base, exons,tx_strand),
        ) |>
        distinct()
    temp[[length(temp) + 1]] <- tx_novel_j
}
tx_novel_junctions <- bind_rows(temp) |>
    left_join(mean_cov, join_by(tx_id == gene, sample_id == sample)) |>
    left_join(selected.genes.and.transcripts, join_by(tx_id == transcript_id))

# Summarize worst ones
# junction_ratio is the ratio of number of unique reads spanning the junction to the mean coverage of the gene
tx_novel_junctions |>
    mutate(junction_ratio = max_num_unique_reads_spanning / mean_cov) |>
    filter(max_num_unique_reads_spanning > 5) |>
    group_by(gene_id, gene_name, tx_id, first_base, last_base, first_base_tx, last_base_tx) |>
    summarize(max_junction_ratio = max(junction_ratio), .groups="drop") |>
    #select(gene_id, gene_name, tx_id, junction_ratio, first_base, last_base, first_base_tx, last_base_tx) |>
    arrange(max_junction_ratio) |>
    write_tsv(output)

#junction_info <- tx_novel_junctions |> left_join(
#    mean_cov |> group_by(gene) |> summarize(min_mean_cov = min(mean_cov)),
#    join_by(tx_id == gene)
#) |> select(tx_id, contig, min_mean_cov, max_num_unique_reads_spanning)
