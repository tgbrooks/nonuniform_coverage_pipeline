library(BSgenome)
library(ensembldb)
library(GenomicRanges)
library(GenomeInfoDb)
library(GenomicAlignments)
library(Rsamtools)
library(tibble)
library(readr)
library(dplyr)

#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#bam.file <- "data/SRX6902444/bam/Aligned.sortedByCoord.out.bam"
#selected_transcripts_file <- "data/mouse_muscle/high_expressed_single_isoform_genes.txt"
#outfile <- "temp.txt"
#sample_id <- "SRX6902444"
#paired <- FALSE
#strandedness <- "reverse"

ensembldb_sqlite <- snakemake@input$ensdb
bam.file <- snakemake@input$bam
selected_transcripts_file <- snakemake@input$transcripts
outfile <- snakemake@output$cov
sample_id <- snakemake@wildcards$sample_id
paired <- snakemake@params$paired
strandedness <- snakemake@params$strandedness


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
        first.plus <- as.logical(strand(GenomicAlignments::first(x)) == "+")
        ifelse(first.plus, start(GenomicAlignments::first(x)), start(GenomicAlignments::last(x)))
    } else {
        start(x)
    }
}
startRight <- function(x) {
    if (inherits(x, "GAlignmentPairs")) {
        first.plus <- as.logical(strand(GenomicAlignments::first(x)) == "+")
        ifelse(first.plus, end(GenomicAlignments::first(x)), end(GenomicAlignments::last(x)))
    } else {
        end(x)
    }
}
endLeft <- function(x) {
    if (inherits(x, "GAlignmentPairs")) {
        first.plus <- as.logical(strand(GenomicAlignments::first(x)) == "+")
        ifelse(first.plus, start(GenomicAlignments::last(x)), start(GenomicAlignments::first(x)))
    } else {
        start(x)
    }
}
endRight <- function(x) {
    if (inherits(x, "GAlignmentPairs")) {
        first.plus <- as.logical(strand(GenomicAlignments::first(x)) == "+")
        ifelse(first.plus, end(GenomicAlignments::last(x)), end(GenomicAlignments::first(x)))
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
gaToReadsOnTx <- function(ga, grl, fco=NULL, connect_reads=TRUE) {
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
        if (connect_reads) {
            # depending on strand of gene:
            # start of left will be the first coordinate on the transcript (+ gene)
            # or start of left will be the last coordinate on the transcript (- gene)
            if (strand %in% c("+", "*")) {
                start <- genomeToTx(startLeft(ga[read.idx]), map)
                end <- genomeToTx(endRight(ga[read.idx]), map)
            } else if (strand == "-") {
                start <- genomeToTx(endRight(ga[read.idx]), map)
                end <- genomeToTx(startLeft(ga[read.idx]), map)
            } else {
            }
            valid <- start < end & !is.na(start) & !is.na(end)
            reads[[length(reads)+1]] <- IRanges(start[valid], end[valid])
        } else {
            # depending on strand of gene:
            # start of left will be the first coordinate on the transcript (+ gene)
            # or start of left will be the last coordinate on the transcript (- gene)
            if (strand %in% c("+", "*")) {
                start1 <- genomeToTx(startLeft(ga[read.idx]), map)
                end1 <- genomeToTx(startRight(ga[read.idx]), map)
                start2 <- genomeToTx(endLeft(ga[read.idx]), map)
                end2 <- genomeToTx(endRight(ga[read.idx]), map)
            } else if (strand == "-") {
                start1 <- genomeToTx(endRight(ga[read.idx]), map)
                end1 <- genomeToTx(endLeft(ga[read.idx]), map)
                start2 <- genomeToTx(startRight(ga[read.idx]), map)
                end2 <- genomeToTx(startLeft(ga[read.idx]), map)
            }
            # Add first (5' most) read
            valid <- start1 < end1 & !is.na(start1) & !is.na(end1)
            reads1 <- IRanges(start1[valid], end1[valid])
            # Add second (3' most) read if not completely redundant
            newstart2 <- pmax(start2, end1+1) # Remove overlap with first read
            valid <- newstart2 < end2 & !is.na(newstart2) & !is.na(end2)
            reads2 <- IRanges(newstart2[valid], end2[valid])
            reads[[i]] <- c(reads1, reads2)
        }
    }
    names(reads) <- names(grl)
    reads
}

getCoverage <- function(gene, bam.file, nonunique=FALSE, connect_reads=TRUE, paired=TRUE, strandedness="unstranded") {
    #stopifnot(is(gene, "GRanges"))

    generange <- range(gene)
    strand(generange) <- "*" # not necessary

    if (nonunique) {
        param <- ScanBamParam(which=generange, what=c("flag","mrnm","mpos"), flag=scanBamFlag())
    } else {
        param <- ScanBamParam(which=generange, what=c("flag","mrnm","mpos"), flag=scanBamFlag(isSecondaryAlignment=FALSE))
    }

    if (strandedness == "forward") {
        strand_mode = 1
    } else if (strandedness == "reverse") {
        strand_mode = 2
    } else {
        strand_mode = 0
        strand(gene) <- "*" # Necessary for findCompatibleOverlaps to function with unstranded reads
    }

    ga <- readGAlignments(bam.file, use.names=TRUE, param=param)
    if (paired) {
        ga <- makeGAlignmentPairs(ga, strandMode=strand_mode)
    } else {
        if (strandedness == "reverse") {
            # Invert the case
            strand(ga) <- case_when(
                decode(strand(ga)) == "+" ~ "-",
                decode(strand(ga)) == "-" ~ "+",
                decode(strand(ga)) == "*" ~ "*",
            )
        } else if(strandedness == "unstranded") {
            strand(ga) <- "*"
        }
    }

    # Check strandedness
    gene_strand <- decode(strand(gene))[[1]] # use strand of first exon - all should be the same
    strand_compatible <- gene_strand == strand(ga)
    ga <- ga[strand_compatible]

    ga <- GenomeInfoDb::keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    fco <- GenomicAlignments::findCompatibleOverlaps(ga, GRangesList(gene))
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco, connect_reads)

    l <- sum(width(gene))
    frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])

    # add 0s out to the length of the gene if no reads go all the way
    c(frag.cov, rep(0, l-length(frag.cov)))
}


res <- list()
for (gene_name in names(ebt.fit)) {
    cov <- getCoverage(
        gene=ebt.fit[[gene_name]],
        bam.file=bam.file,
        nonunique=FALSE,
        connect_reads=TRUE,
        paired=paired,
        strandedness=strandedness
    )
    cov_nonuniq <- getCoverage(
        gene=ebt.fit[[gene_name]],
        bam.file=bam.file,
        nonunique=TRUE,
        connect_reads=TRUE,
        paired=paired,
        strandedness=strandedness
    )
    cov_disconnected <- getCoverage(
        gene=ebt.fit[[gene_name]],
        bam.file=bam.file,
        nonunique=FALSE,
        connect_reads=FALSE,
        paired=paired,
        strandedness=strandedness
    )
    res[[gene_name]] <- tibble(
        sample = sample_id,
        gene = gene_name,
        pos = 1:length(cov),
        cov = decode(cov),
        cov_nonunique = decode(cov_nonuniq),
        cov_disconnected = decode(cov_disconnected),
    )
}
all <- bind_rows(res)
all |> write_tsv(outfile)
