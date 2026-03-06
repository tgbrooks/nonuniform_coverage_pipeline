library(ensembldb)
library(readr)
library(dplyr)
library(ggplot2)

sample_ids <- snakemake@params$sample_ids
ensembldb_sqlite <- snakemake@input$ensdb
outdir <- snakemake@output$outdir
cov_files <- snakemake@input$cov_tables
transcript_file <- snakemake@input$transcripts
sample_info_file <- snakemake@input$sample_info
tissue <- snakemake@wildcards$tissue

# TESTING VALUES ###
#sample_ids <- c("SRX4080514", "SRX4080520", "SRX13396186", "SRX13396189") 
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#transcript_file <- "data/liver/high_expressed_single_isoform_genes.txt"
#outdir <- "results/alpine_fit_plots/temp/"
#cov_files <- paste0("results/alpine_fits/liver/", sample_ids, ".coverage_table.txt")
#sample_info_file <- "results/liver.sample_info.txt"
#tissue <- "liver"
###########

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
if (tissue == "liver") {
    select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973")
} else if (tissue == "UHR") {
    select_genes <- c("ENST00000380680", "ENST00000604000", "ENST00000426077")
} else if (tissue == "UHR_degraded") {
    select_genes <- c("ENST00000053468", "ENST00000253063", "ENST00000272233")
} else {
    select_genes <- c()
}

if (length(select_genes) > 0) {
    print(select_genes)
    select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
        left_join(sample_info, by=join_by(sample == ID)) %>%
        group_by(sample, gene) %>%
        mutate(rel_read_depth = actual / max(actual)) %>%
        arrange(study) %>%
        ungroup() %>%
        left_join(high_exp_genes, by=join_by(gene == transcript_id)) %>%
        mutate(full_gene = paste0(gene_name, "\n", gene))
    exon_info <- lapply(
            select_genes,
            function(gene) { get_exon_cds_info(gene) |> mutate(gene = gene) }
        ) |>
        bind_rows() |>
        left_join(high_exp_genes, join_by(gene == transcript_id)) |>
        mutate(full_gene = paste0(gene_name, "\n", gene)) |>
        cross_join(
            # need one copy per study
            sample_info |>
                filter(study %in% select_cov$study) |>
                select(study) |>
                distinct()
        )
    ggplot(
            data = select_cov,
            aes(x=pos / 1000, y=rel_read_depth)
        ) +
        facet_grid(
            rows=vars(study),
            cols=vars(full_gene),
            scales = "free",
        ) +
        geom_path(
            aes(color = sample_num),
        ) +
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
            color = "Sample\nnumber"
        ) +
        scale_fill_manual(values=c("#888", "#444")) + 
        theme(
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    ggsave(
        paste(outdir, "selected_genes.replicates.png", sep="/"),
        width = 7,
        height = 7,
    )
}

if (("rin_score" %in% colnames(sample_info)) & (any(!is.na(sample_info$rin_score)))) {
    # BY RIN SCORE
    select_cov <- cov_table[cov_table$gene %in% select_genes,] %>%
        left_join(sample_info, by=join_by(sample == ID)) %>%
        group_by(sample, gene) %>%
        mutate(rel_read_depth = actual / max(actual)) %>%
        ungroup() %>%
        left_join(high_exp_genes, by=join_by(gene == transcript_id)) %>%
        mutate(full_gene = paste0(gene_name, "\n", gene))
    ggplot(
            data = select_cov,
            aes(x=pos / 1000, y=rel_read_depth, group=sample)
        ) +
        facet_grid(
            cols = vars(full_gene),
            scales = "free",
        ) +
        geom_path(
            aes(color = rin_score),
        ) +
        scale_color_viridis_c(option = "viridis", limits=c(0,10)) +
        labs(
            x = "Position (kb)",
            y = "Normalized read depth",
            color = "RIN score"
        ) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    ggsave(
        paste(outdir, "selected_genes.by_rin_score.png", sep="/"),
        width = 7,
        height = 2.5,
    )
}

# Make the plots
for (gene in unique(cov_table$gene)) {
  print(paste("Plotting", gene))
  gene_cov <- cov_table[cov_table$gene == gene,] %>%
        left_join(sample_info, by=join_by(sample == ID))
  print(gene_cov |> colnames())
  ggplot(
          data = gene_cov,
          aes(x=pos / 1000, y=actual)
      ) +
      facet_grid(
          cols=vars(sample_num),
          rows=vars(study),
          scales = "free",
      ) +
      geom_path(
          color = "black",
      ) +
      geom_path(
          aes(y = predicted),
          color = "red"
      ) +
      labs(
          x = "Position (kb)",
          y = "Read depth",
          title = "Alpine model fits"
      )
  ggsave(
      paste(outdir, paste("trace.", gene, ".png", sep=""), sep="/"),
      width = 7,
      height = 7,
  )
}
