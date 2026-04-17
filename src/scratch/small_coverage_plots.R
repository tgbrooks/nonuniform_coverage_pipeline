library(ensembldb)
library(tidyverse)

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


### MOUSE STUDIES
txdb <- EnsDb("data/Mus_musculus.GRCm38.102.gtf.sqlite")
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
ebt0 <- exonsBy(txdb, by="tx")


## LIVER
sample_ids = c(
    "SRX4080514",
    "SRX4080520",
    #"SRX6902444",
    #"SRX6902450",
    #"SRX3304756",
    #"SRX3304764",
    "SRX4393368",
    "SRX4393369",
    "SRX16386863",
    "SRX16386864",
    "SRX14468350",
    "SRX14468347",
    "SRX13396189",
    "SRX13396186",
    "SRX11694510",
    "SRX11694499"
)
tissue <- "liver"
high_exp_genes <- read_tsv("data/liver/high_expressed_single_isoform_genes.txt")

cov_tables_files = paste0("data/", sample_ids, "/transcript_coverage.txt")

cov_tables <- list()
for (file in cov_tables_files) {
  message("Reading ", file)
  cov_tables[[length(cov_tables)+1]] <- read_tsv(file, show_col_types=FALSE) |>
      select(
             sample_id = sample,
             transcript_id = gene,
             pos = pos,
             cov = cov,
    )
}
cov_table <- bind_rows(cov_tables)

sample_info <- read_tsv(paste0("results/", tissue, ".sample_info.txt")) %>%
    group_by(study) %>%
    mutate(sample_num = row_number()) %>%
    ungroup() %>%
    mutate(sample_num = as.factor(sample_num))

# Plot for the paper
#select_genes <- c("ENSMUST00000023559", "ENSMUST00000028995", "ENSMUST00000047973") # old paper v1 values
select_genes <- c("ENSMUST00000025218")#, "ENSMUST00000025356", "ENSMUST00000027144") # for paper v2
selected_studies <- c("GSE117134", "PRJNA816471")
select_cov <- cov_table %>%
    filter(transcript_id %in% select_genes) %>%
    left_join(sample_info, by=join_by(sample_id == ID)) %>%
    filter(study %in% selected_studies) %>%
    group_by(sample_id, transcript_id) %>%
    mutate(cov_norm = cov / max(cov)) %>%
    left_join(high_exp_genes, by=join_by(transcript_id == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", transcript_id)) %>%
    arrange(study) %>%
    ungroup()
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
        aes(x=pos / 1000, y=cov_norm)
    ) +
    facet_grid(
        cols = vars(full_gene),
        scales = "free",
    ) +
    geom_path(
        aes(color = study, group=sample_id, linewidth=sample_num),
    ) +
    scale_linewidth_manual(values=c("1"=1.0, "2"=0.5), name='Sample') +
    #scale_alpha_discrete(range=c(0.5,1.0), name="sample")+
    # THE EXON ANNOTATION LAYERS
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
    scale_fill_manual(values=c("#888", "#444")) + 
    labs(
        x = "Position (kb)",
        y = "Normalized coverage",
        color = "Study"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "small.liver.selected_genes.replicates.png", sep="/"),
    width = 5,
    height = 3,
)

#### HUMAN
txdb <- EnsDb("data/GRCh38.ensemblv109.gtf.sqlite")
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
ebt0 <- exonsBy(txdb, by="tx")

## UHR
sample_ids = c(
    "SRX302130",
    "SRX302146",
    "SRX302162",
    "SRX302178",
    "SRX302514",
    "SRX302529",
    "SRX302544",
    "SRX302559",
    "SRX302874",
    "SRX302890",
    "SRX302906",
    "SRX302922",
    "SRX302574",
    "SRX302194",
    "SRX302938"
)
tissue <- "UHR"

high_exp_genes <- read_tsv("data/UHR/high_expressed_single_isoform_genes.txt")

cov_tables_files = paste0("data/", sample_ids, "/transcript_coverage.txt")
cov_tables <- list()
for (file in cov_tables_files) {
  message("Reading ", file)
  cov_tables[[length(cov_tables)+1]] <- read_tsv(file, show_col_types=FALSE) |>
      select(
             sample_id = sample,
             transcript_id = gene,
             pos = pos,
             cov = cov,
    )
}
cov_table <- bind_rows(cov_tables)

sample_info <- read_tsv(paste0("results/", tissue, ".sample_info.txt")) %>%
    mutate(
           site = str_match(study, "_([A-Z]+)")[,2]
        )

# Plot for the paper
select_genes <- c("ENST00000269593") #choice for v2
selected_sites <- c("BGI", "MAY")
select_cov <- cov_table %>%
    filter(transcript_id %in% select_genes) %>%
    left_join(sample_info, by=join_by(sample_id == ID)) %>%
    filter(site %in% selected_sites) %>%
    group_by(sample_id, transcript_id) %>%
    mutate(cov_norm = cov / max(cov)) %>%
    left_join(high_exp_genes, by=join_by(transcript_id == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", transcript_id)) %>%
    arrange(site) %>%
    ungroup()
exon_info <- lapply(
        select_genes,
        function(gene) { get_exon_cds_info(gene) |> mutate(gene = gene) }
    ) |>
    bind_rows() |>
    left_join(high_exp_genes, join_by(gene == transcript_id)) |>
    mutate(full_gene = paste0(gene_name, "\n", gene)) |>
    cross_join(
        # need one copy per site
        sample_info |>
            filter(site %in% select_cov$site) |>
            select(site) |>
            distinct()
    )
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=cov_norm)
    ) +
    facet_grid(
        cols = vars(full_gene),
        scales = "free",
    ) +
    geom_path(
        aes(color = site, group=sample_id, linewidth=library_id),
    ) +
    scale_linewidth_continuous(range=c(0.3,1.0), name='Sample') +
    # THE EXON ANNOTATION LAYERS
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
    scale_fill_manual(values=c("#888", "#444")) + 
    labs(
        x = "Position (kb)",
        y = "Normalized coverage",
        color = "Site"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "small.UHR.selected_genes.replicates.png", sep="/"),
    width = 5,
    height = 3,
)


### UHR_degraded
tissue <- "UHR_degraded"
high_exp_genes <- read_tsv("data/UHR_degraded/high_expressed_single_isoform_genes.txt")
sample_info <- read_tsv(paste0("results/", tissue, ".sample_info.txt"))
sample_ids <- sample_info$ID
cov_tables_files = paste0("data/", sample_ids, "/transcript_coverage.txt")
cov_tables <- list()
for (file in cov_tables_files) {
  message("Reading ", file)
  cov_tables[[length(cov_tables)+1]] <- read_tsv(file, show_col_types=FALSE) |>
      select(
             sample_id = sample,
             transcript_id = gene,
             pos = pos,
             cov = cov,
    )
}
cov_table <- bind_rows(cov_tables)

select_genes <- c("ENST00000272233") #choice for v2

# BY RIN SCORE
select_cov <- cov_table %>%
    filter(transcript_id %in% select_genes) %>%
    left_join(sample_info, by=join_by(sample_id == ID)) %>%
    group_by(sample_id, transcript_id) %>%
    mutate(rel_read_depth = cov / max(cov)) %>%
    ungroup() %>%
    left_join(high_exp_genes, by=join_by(transcript_id == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", transcript_id))
exon_info <- lapply(
        select_genes,
        function(gene) { get_exon_cds_info(gene) |> mutate(gene = gene) }
    ) |>
    bind_rows() |>
    left_join(high_exp_genes, join_by(gene == transcript_id)) |>
    mutate(full_gene = paste0(gene_name, "\n", gene))
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=rel_read_depth, group=sample_id)
    ) +
    facet_grid(
        cols = vars(full_gene),
        scales = "free",
    ) +
    geom_path(
        aes(color = rin_score),
    ) +
    scale_color_viridis_c(option = "viridis", limits=c(0,10)) +
    # THE EXON ANNOTATION LAYERS
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
    scale_fill_manual(values=c("#888", "#444")) + 
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
    paste("results", "scratch", "small.UHR_degraded.selected_genes.png", sep="/"),
    width = 5,
    height = 3,
)
