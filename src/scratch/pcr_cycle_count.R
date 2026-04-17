library(ensembldb)
library(tidyverse)
### GVS for UHR_PCR dataset
uhr <- read_tsv("results/UHR_PCR.sample_info.txt")
temp <- list()
for (sample_id in uhr$ID) {
    temp[[length(temp) + 1]] <- read_tsv(paste0("data/", sample_id, "/corrected_global_valley_score.txt"))
}

gvs <- bind_rows(temp) |>
    left_join(
        uhr,
        join_by(sample == ID),
    )

ggplot(gvs, aes(x=PCR_cycle_count, y = global_valley_score, color=input_material_ng))  +
    scale_color_viridis_c(trans="log10") +
    labs(
         x = "PCR cycle number",
         y = "Global valley score",
         color = "Input material (ng)"
    ) +
    #facet_wrap(vars(input_material_ng)) +
    geom_point()
ggsave("results/scratch/GVS.by_PCR_cycle_count.png", width=5, height=5)
PCR_cycle_pvalue <- summary(lm(global_valley_score ~ PCR_cycle_count, data=gvs))$coefficients[2,4]
message(paste("GVS ~ PCR cycle count pvalue", PCR_cycle_pvalue))


### Coverage plot by PCR cycle number
ensembldb_sqlite <- "data/GRCh38.ensemblv109.gtf.sqlite"

sample_ids = uhr$ID

# Load coverage
temp <- list()
for (sample_id in sample_ids) {
    temp[[length(temp)+1]] <- read_tsv(paste0("data/", sample_id, "/transcript_coverage.txt"), show_col_types=FALSE)
}
cov <- bind_rows(temp)


# Load annotations
txdb <- EnsDb(ensembldb_sqlite)
genedf <- genes(txdb, return.type="DataFrame")
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

ebt0 <- exonsBy(txdb, by="tx")

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
        length = lengths(tx),
        exon_id = tx$exon_id,
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
        ) |>
        mutate(group = paste(type, parity))
    return(exons2)
}

# Plot for the paper
#select_genes <- c("ENST00000307630", "ENST00000340857", "ENST00000272233")
#select_genes <- c("ENST00000361624", "ENST00000217407", "ENST00000307630")
select_genes <- c("ENST00000361899", "ENST00000361381",  "ENST00000361624")

sample_info <- uhr %>%
    mutate(
           sample_id = ID
           )
print(select_genes)
select_cov <- cov[cov$gene %in% select_genes,] %>%
    left_join(sample_info, by=join_by(sample == sample_id)) %>%
    group_by(sample, gene) %>%
    mutate(rel_read_depth = cov / max(cov)) %>%
    ungroup() %>%
    left_join(txdf %>% as_tibble(), by=join_by(gene == tx_id)) %>%
    left_join(genedf %>% as_tibble(), by=join_by(gene_id==gene_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", gene))
exon_info <- lapply(
        select_genes,
        function(gene) { get_exon_cds_info(gene) |> mutate(gene = gene) }
    ) |>
    bind_rows() |>
    left_join(txdf %>% as_tibble(), by=join_by(gene == tx_id)) %>%
    left_join(genedf %>% as_tibble(), by=join_by(gene_id==gene_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", gene))

ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=rel_read_depth)
    ) +
    facet_grid(
        cols=vars(full_gene),
        rows=vars(input_material_ng),
        scales = "free",
        labeller = labeller(
            input_material_ng = function(x) paste(x, "ng")
        )
    ) +
    geom_path(
        aes(color = PCR_cycle_count, group=sample),
    ) +
    scale_color_viridis_c(option="A") +
    # the exon annotation layer
    geom_rect(
        aes(
            xmin = exon_start/1000,
            xmax = exon_end/1000,
            ymin=case_match(type, "utr"~-0.075, "cds"~-0.1),
            ymax=case_match(type, "utr"~-0.025, "cds"~0.0),
            fill=as.factor(parity)),
        data = exon_info,
        show.legend=FALSE,
        inherit.aes=FALSE,
    ) + 
    labs(
        x = "Position (kb)",
        y = "Normalized read depth",
        color = "PCR cycle number"
    ) +
    scale_fill_manual(values=c("#888", "#444")) + 
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    "results/scratch/UHR_PCR_seq_coverage.png",
    width = 7,
    height = 10,
)
