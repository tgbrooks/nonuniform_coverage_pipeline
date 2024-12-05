library(BSgenome)
library(ensembldb)
library(GenomicRanges)
library(GenomeInfoDb)
library(GenomicAlignments)
library(Rsamtools)
library(tibble)
library(readr)
library(dplyr)

#ensembldb_sqlite <- "data/GRCh38.ensemblv109.gtf.sqlite"
#bam.file <- "data/SRX554541/bam/Aligned.sortedByCoord.out.bam"
#selected_transcripts_file <- "data/UHR/high_expressed_single_isoform_genes.txt"
#outfile <- "temp.txt"
#sample_id <- "SRX554541"

ensembldb_sqlite <- snakemake@input$ensdb
bam.file <- snakemake@input$bam
selected_transcripts_file <- snakemake@input$transcripts
outfile <- snakemake@output$cov
sample_id <- snakemake@wildcards$sample_id

txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

selected.genes.and.transcripts <- read.csv(selected_transcripts_file, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]


# Adapted from mikelove/alpine
startLeft <- function(x) {
    if (inherits(x, "GAlignmentPairs")) {
        first.plus <- as.logical(strand(first(x)) == "+")
        ifelse(first.plus, start(first(x)), start(last(x)))
    } else {
        start(x)
    }
}
endRight <- function(x) {
    if (inherits(x, "GAlignmentPairs")) {
        first.plus <- as.logical(strand(first(x)) == "+")
        ifelse(first.plus, end(last(x)), end(first(x)))
    } else {
        end(x)
    }
}
mapTxToGenome <- function(exons) {
    strand <- as.character(strand(exons)[1])
    stopifnot(all(exons$exon_rank == seq_along(exons)))

    # Hack to replicate `rev` from old S4Vectors:::fancy_mseq
    if(strand == "-"){
        froms <- start(exons) + width(exons) - 1L
        bys <- -1L

    }else{
        froms <- start(exons)
        bys <- 1L

    }

    bases <- sequence(width(exons), from = froms,
                      by = bys)

    data.frame(tx=seq_along(bases),
               genome=bases,
               exon_rank=rep(exons$exon_rank, width(exons)))
}
genomeToTx <- function(genome, map) map$tx[match(genome, map$genome)]
txToGenome <- function(tx, map) map$genome[match(tx, map$tx)]
txToExon <- function(tx, map) map$exon_rank[match(tx, map$tx)]
gaToReadsOnTx <- function(ga, grl, fco=NULL) {
    reads <- list()
    for (i in seq_along(grl)) {
        exons <- grl[[i]]
        strand <- as.character(strand(exons)[1])
        read.idx <- if (is.null(fco)) {
            seq_along(ga)
        } else {
            queryHits(fco)[subjectHits(fco) == i]
        }
        map <- mapTxToGenome(exons)
        # depending on strand of gene:
        # start of left will be the first coordinate on the transcript (+ gene)
        # or start of left will be the last coordinate on the transcript (- gene)
        if (strand == "+") {
            start <- genomeToTx(startLeft(ga[read.idx]), map)
            end <- genomeToTx(endRight(ga[read.idx]), map)
        } else if (strand == "-") {
            start <- genomeToTx(endRight(ga[read.idx]), map)
            end <- genomeToTx(startLeft(ga[read.idx]), map)
        }
        valid <- start < end & !is.na(start) & !is.na(end)
        reads[[i]] <- IRanges(start[valid], end[valid])
    }
    names(reads) <- names(grl)
    reads
}

getCoverage <- function(gene, bam.file) {
    stopifnot(is(gene, "GRanges"))

    generange <- range(gene)
    strand(generange) <- "*" # not necessary

    param <- ScanBamParam(which=generange, what=c("flag","mrnm","mpos"), flag=scanBamFlag(isSecondaryAlignment=FALSE))
    # TODO: this assume single-end reads
    ga <- readGAlignments(bam.file, use.names=TRUE, param=param)

    ga <- GenomeInfoDb::keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    fco <- GenomicAlignments::findCompatibleOverlaps(ga, GRangesList(gene))
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)

    l <- sum(width(gene))
    frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])

    # add 0s out to the length of the gene if no reads go all the way
    c(frag.cov, rep(0, l-length(frag.cov)))
}


res <- list()
for (gene in names(ebt.fit)) {
    cov <- getCoverage(
        gene=ebt.fit[[gene]],
        bam.file=bam.file
    )
    res[[gene]] <- tibble(
        sample = sample_id,
        gene = gene,
        pos = 1:length(cov),
        cov = decode(cov),
    )
}
all <- bind_rows(res)
all |> write_tsv(outfile)
