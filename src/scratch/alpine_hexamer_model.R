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

get_exon_cds_info <- function(tx_id, ebt0) {
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

#special_model_types <- c("hexamer", "pos_no_genelen", "pos_only", "GC", "hexamer_no_pos", "baseline")
special_model_types <- c("hexamer", "GC", "pos_only")
sample_id <- "SRX11694510"
model_cov <- lapply(
       special_model_types,
       function(model_type) {
    read_tsv(paste0("results/alpine_fits", model_type, "/liver/", sample_id, ".coverage_table.txt")) |>
        select(
            sample_id = sample,
            transcript_id = gene,
            pos = pos,
            cov = predicted,
        ) |> mutate(
            type = model_type 
        )
       }
       )
cov <- bind_rows(
    read_tsv(paste0("results/alpine_fits/liver/", sample_id, ".coverage_table.txt")) |>
        select(
            sample_id = sample,
            transcript_id = gene,
            pos = pos,
            cov = actual,
        ) |> mutate(
            type = 'actual'
        ),
    model_cov,
    read_tsv(paste0("results/alpine_fits/liver/", sample_id, ".coverage_table.txt")) |>
        select(
            sample_id = sample,
            transcript_id = gene,
            pos = pos,
            cov = predicted,
        ) |> mutate(
            type = 'full'
        ),
)


############# LIVER DATA ###################
txdb <- EnsDb("data/GRCh38.ensemblv109.gtf.sqlite")
hg_txdf <- transcripts(txdb, return.type="DataFrame")

txdb <- EnsDb("data/Mus_musculus.GRCm38.102.gtf.sqlite")
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
ebt0 <- exonsBy(txdb, by="tx")

high_exp_genes <- read_tsv("data/liver/high_expressed_single_isoform_genes.txt")
### PLOT COVERAGE VS ALPINE FIT
#select_genes <- c("ENSMUST00000025218", "ENSMUST00000025356", "ENSMUST00000027144") # for paper v2
#select_genes <- c("ENSMUST00000025218", "ENSMUST00000021001", "ENSMUST00000023952")
select_genes <- c("ENSMUST00000044355", "ENSMUST00000030538", "ENSMUST00000082059") # not in training set
select_cov <- cov %>%
    filter(transcript_id %in% select_genes) %>%
    group_by(sample_id, transcript_id) %>%
    left_join(high_exp_genes, by=join_by(transcript_id == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", transcript_id)) %>%
    ungroup()
exon_info <- lapply(
        select_genes,
        function(gene) { get_exon_cds_info(gene, ebt0) |> mutate(gene = gene) }
    ) |>
    bind_rows() |>
    left_join(high_exp_genes, join_by(gene == transcript_id)) |>
    mutate(full_gene = paste0(gene_name, "\n", gene))
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=cov, color=type)
    ) +
    facet_grid(
        #rows=vars(study),
        cols=vars(full_gene),
        scales = "free",
    ) +
    geom_path() +
    scale_color_manual(values=c("actual"="black", "hexamer"="orange", "full"="red", "pos_no_genelen"="blue", "pos_only"="navy", "hexamer_no_pos"="yellow", "GC"="green"), name="") +
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
        y = "Depth",
        color = "Study"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "alpine_hexamer_model.png", sep="/"),
    width = 9,
    height = 5,
)

############ IVT-SEQ DATA ##########################
sample_id <- "SRX341390"
model_cov <- lapply(
       special_model_types,
       function(model_type) {
    read_tsv(paste0("results/alpine_fits", model_type, "/IVT_seq/", sample_id, ".coverage_table.txt")) |>
        select(
            sample_id = sample,
            transcript_id = gene,
            pos = pos,
            cov = predicted,
        ) |> mutate(
            type = model_type 
        )
       }
       )
cov <- bind_rows(
    read_tsv(paste0("results/alpine_fitshexamer/IVT_seq/", sample_id, ".coverage_table.txt")) |>
        select(
            sample_id = sample,
            transcript_id = gene,
            pos = pos,
            cov = actual,
        ) |> mutate(
            type = 'actual'
        ),
    model_cov,
#    read_tsv(paste0("results/alpine_fits/IVT_seq/", sample_id, ".coverage_table.txt")) |>
#        select(
#            sample_id = sample,
#            transcript_id = gene,
#            pos = pos,
#            cov = predicted,
#        ) |> mutate(
#            type = 'full'
#        ),
)
txdb <- EnsDb("data/IVT_seq.gtf.sqlite")
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
ebt0 <- exonsBy(txdb, by="tx")

high_exp_genes <- read_tsv("data/IVT_seq/high_expressed_single_isoform_genes.txt")
select_genes <- c("BC011380", "BC013581", "BC016145")# first one was used in Alpine paper, all from test set
gene_mean <- cov %>% filter(type == "actual") %>% group_by(transcript_id) %>% summarize(mean_cov = mean(cov))
select_cov <- cov %>%
    filter(transcript_id %in% select_genes) %>%
    left_join(gene_mean, "transcript_id") %>%
    group_by(sample_id, transcript_id) %>%
    left_join(high_exp_genes, by=join_by(transcript_id == transcript_id)) %>%
    mutate(full_gene = paste0(gene_name, "\n", transcript_id), cov_norm = cov / mean_cov) %>%
    ungroup()
#exon_info <- lapply(
#        select_genes,
#        function(gene) { get_exon_cds_info(gene, ebt0) |> mutate(gene = gene) }
#    ) |>
#    bind_rows() |>
#    left_join(high_exp_genes, join_by(gene == transcript_id)) |>
#    mutate(full_gene = paste0(gene_name, "\n", gene))
ggplot(
        data = select_cov,
        aes(x=pos / 1000, y=cov_norm, color=type)
    ) +
    facet_grid(
        #rows=vars(study),
        cols=vars(full_gene),
        scales = "free",
    ) +
    geom_path() +
    scale_color_manual(values=c("actual"="black", "hexamer"="orange", "full"="red", "pos_no_genelen"="blue", "pos_only"="navy", "hexamer_no_pos"="yellow", "GC"="green", "baseline"="grey"), name="") +
    #scale_alpha_discrete(range=c(0.5,1.0), name="sample")+
    ## THE EXON ANNOTATION LAYERS
    #geom_rect(
    #    aes(
    #        xmin = exon_start/1000,
    #        xmax = exon_end/1000,
    #        ymin=case_match(type, "utr"~-0.075, "cds"~-0.1),
    #        ymax=case_match(type, "utr"~-0.025, "cds"~0.0),
    #        fill=as.factor(parity)),
    #    data = exon_info,
    #    show.legend=FALSE,
    #    inherit.aes=FALSE,
    #) +
    #scale_fill_manual(values=c("#888", "#444")) +
    labs(
        x = "Position (kb)",
        y = "Depth",
        color = "Study"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(
    paste("results", "scratch", "alpine_hexamer_model.IVT_seq.png", sep="/"),
    width = 9,
    height = 2.5,
)

