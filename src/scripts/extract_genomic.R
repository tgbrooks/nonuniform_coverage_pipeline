# Extract genomic sequence of the selected transcripts for ease of use
library(ensembldb)
library(alpine)
library(readr)
library(dplyr)
library(tibble)
library(GenomicRanges)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "NCBI"
seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "NCBI"

GRCm38_txdb <- EnsDb("data/Mus_musculus.GRCm38.102.gtf.sqlite")
GRCh38_txdb <- EnsDb("data/GRCh38.ensemblv109.gtf.sqlite")

liver_selected <- read_tsv("data/liver/high_expressed_single_isoform_genes.txt")
UHR_selected <- read_tsv("data/UHR/high_expressed_single_isoform_genes.txt")

get_seqs <- function(selected, genome, txdb) {
    selected.txs = selected$transcript_id
    ebt0 <- exonsBy(txdb, by="tx")
    ebt.fit <- ebt0[selected.txs]
    tx.dna <- list()
    for (txid in names(ebt.fit)) {
        exons <- ebt.fit[[txid]]
        map <- alpine:::mapTxToGenome(exons)
        exon.dna <- getSeq(genome, exons)
        tx.dna[[txid]] <- unlist(exon.dna) |> as.character()
    }

    tibble(
        transcript_id = names(tx.dna),
        sequence = unlist(tx.dna),
    )
}

liver_tx_dna <- get_seqs(liver_selected, BSgenome.Mmusculus.UCSC.mm10, GRCm38_txdb)
UHR_tx_dna <- get_seqs(UHR_selected, BSgenome.Hsapiens.UCSC.hg38, GRCh38_txdb)


all_dna <- bind_rows(liver_tx_dna, UHR_tx_dna)

write_tsv(all_dna, "results/selected_genomic_sequences.txt")
