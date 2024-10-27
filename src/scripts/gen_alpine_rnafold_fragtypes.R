#options(error = function() traceback(10))

library(BSgenome)
library(ensembldb)
library(stringr)
library(r2r)

#sample_id <- "SRX16386863"
#readlength <- 151
#modelfile <- "results/alpine_rnafold/SRX16386863/model.rda"
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#transcript_file <- "data/liver/high_expressed_single_isoform_genes.txt"

RNA_FOLD_RESOLUTION <- 25

sample_id <- snakemake@wildcards$sample_id
readlength <- snakemake@params$readlength 
modelfile <- snakemake@output$modelfile
ensembldb_sqlite <- snakemake@input$ensdb
transcript_file <- snakemake@input$transcripts

txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

selected.genes.and.transcripts <- read.csv(transcript_file, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]

library(GenomicRanges)
library(alpine)

gene.lengths <- sum(width(ebt.fit))

minsize <- 70
maxsize <- 600
gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

# Load the BSgenome for our species
library(snakemake@params$BSgenome)
if (snakemake@params$BSgenome == "BSgenome.Hsapiens.UCSC.hg38") {
    my_genome <- BSgenome.Hsapiens.UCSC.hg38
} else {
    my_genome <- BSgenome.Mmusculus.UCSC.mm10
}

# Helper function to run RNAfold and extract the entropy
RNAfold_executable <- "/home/thobr/ViennaRNA/bin/RNAfold"
RNAfold <- function(dna) {
  res <- system2(
      RNAfold_executable,
      args = c("--noPS", "--temp", "42"), # NOTE: 42 degrees Celsius is the RT temperature for TruSeq
      stdout = TRUE,
      input = dna,
  )
  mfe <- as.double(str_match(res[[2]], "[()\\.]+ \\( *([-0-9\\.]+)\\)")[,2])
  return(mfe)
}

# Our modified buildFragtypes to include MFE/nt
buildFragtypes <- function (exons, genome, readlength, minsize, maxsize, gc = TRUE,
    gc.str = TRUE, vlmm = TRUE, RNAfold=TRUE)
{
    stopifnot(is(exons, "GRanges"))
    stopifnot(is(genome, "BSgenome"))
    stopifnot(is.numeric(minsize) & is.numeric(maxsize) & is.numeric(readlength))
    stopifnot(sum(width(exons)) >= maxsize)
    stopifnot(all(c("exon_rank", "exon_id") %in% names(mcols(exons))))
    stopifnot(!any(strand(exons) == "*"))
    npre <- 8
    npost <- 12
    map <- alpine:::mapTxToGenome(exons)
    l <- nrow(map)
    strand <- as.character(strand(exons)[1])
    start <- rep(seq_len(l - minsize + 1), each = maxsize - minsize +
        1)
    end <- as.integer(start + minsize:maxsize - 1)
    mid <- as.integer(0.5 * (start + end))
    relpos <- mid/l
    fraglen <- as.integer(end - start + 1)
    id <- IRanges(start, end)
    fragtypes <- DataFrame(start = start, end = end, relpos = relpos,
        fraglen = fraglen, id = id)
    fragtypes <- fragtypes[fragtypes$end <= l, , drop = FALSE]
    exon.dna <- getSeq(genome, exons)
    tx.dna <- unlist(exon.dna)
    if (vlmm) {
        fragtypes$fivep.test <- fragtypes$start - npre >= 1
        fragtypes$fivep <- as(Views(tx.dna, fragtypes$start -
            ifelse(fragtypes$fivep.test, npre, 0), fragtypes$start +
            npost), "DNAStringSet")
        fragtypes$threep.test <- fragtypes$end + npre <= length(tx.dna)
        fragtypes$threep <- as(Views(tx.dna, fragtypes$end -
            npost, fragtypes$end + ifelse(fragtypes$threep.test,
            npre, 0), ), "DNAStringSet")
        fragtypes$threep <- reverseComplement(fragtypes$threep)
    }
    if (gc) {
        fragrange <- minsize:maxsize
        gc.vecs <- lapply(fragrange, function(i) {
            letterFrequencyInSlidingView(tx.dna, view.width = i,
                letters = "CG", as.prob = TRUE)
        })
        fragtypes <- fragtypes[order(fragtypes$fraglen), , drop = FALSE]
        fragtypes$gc <- do.call(c, gc.vecs)
        fragtypes <- fragtypes[order(fragtypes$start), , drop = FALSE]
    }
    if (gc.str) {
        gc.40 <- as.numeric(letterFrequencyInSlidingView(tx.dna,
            40, letters = "CG", as.prob = TRUE))
        max.gc.40 <- max(Views(gc.40, fragtypes$start, fragtypes$end -
            40 + 1))
        gc.20 <- as.numeric(letterFrequencyInSlidingView(tx.dna,
            20, letters = "CG", as.prob = TRUE))
        max.gc.20 <- max(Views(gc.20, fragtypes$start, fragtypes$end -
            20 + 1))
        fragtypes$GC40.90 <- as.numeric(max.gc.40 >= 36/40)
        fragtypes$GC40.80 <- as.numeric(max.gc.40 >= 32/40)
        fragtypes$GC20.90 <- as.numeric(max.gc.20 >= 18/20)
        fragtypes$GC20.80 <- as.numeric(max.gc.20 >= 16/20)
    }
    if (RNAfold) {
      mfe_store <- hashmap() # Memoized RNAfold MFE values
      mfe_per_nt <- apply(
        fragtypes,
        1,
        function (frag) {
          # Round to the specified resolution
          round_start <- round(frag$start/RNA_FOLD_RESOLUTION)*RNA_FOLD_RESOLUTION
          round_end <- min(round(frag$end/RNA_FOLD_RESOLUTION)*RNA_FOLD_RESOLUTION, length(tx.dna))
          if (is.null(mfe_store[[c(round_start, round_end)]])) {
            # Compute value
            print(paste("Computing for", round_start, round_end))
            frag.dna <- as.character(tx.dna[round_start:round_end])
            mfe_per_nt <- RNAfold(frag.dna) / (round_end - round_start + 1)
            mfe_store[[c(round_start, round_end)]] <- mfe_per_nt
          } else {
            mfe_per_nt <- mfe_store[[c(round_start, round_end)]] # use memoized value
          }
          return(mfe_per_nt)
        }
      )
      fragtypes$MFE_per_nt <- mfe_per_nt
    }
    fragtypes$gstart <- alpine:::txToGenome(fragtypes$start, map)
    fragtypes$gend <- alpine:::txToGenome(fragtypes$end, map)
    fragtypes$gread1end <- alpine:::txToGenome(fragtypes$start + readlength -
        1, map)
    fragtypes$gread2start <- alpine:::txToGenome(fragtypes$end - readlength +
        1, map)
    fragtypes
}

fragtypes <- lapply(
    gene.names, function(gene.name) {
      buildFragtypes(exons=ebt.fit[[gene.name]],
                     genome=my_genome,
                     readlength=readlength,
                     minsize=minsize,
                     maxsize=maxsize,
                     gc.str=TRUE)
    }
)

# Add genelen to the data so that we can model based off it(Alpine provides only fraglen)
# We consider genelength as important particularly for polyA selection, which has a large
# 3' bias in long genes but that effect dies off for short genes
for(gene in names(fragtypes)) {
    fragtypes[[gene]]$genelen <- max(fragtypes[[gene]]$end) - min(fragtypes[[gene]]$start)
}


saveRDS(fragtypes, snakemake@output$fragtypes)
