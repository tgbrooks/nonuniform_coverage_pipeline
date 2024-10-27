#options(error = function() traceback(10))

# Load the BSgenome for our species
library(snakemake@params$BSgenome, character.only=TRUE)
if (snakemake@params$BSgenome == "BSgenome.Hsapiens.UCSC.hg38") {
    my_genome <- BSgenome.Hsapiens.UCSC.hg38
} else {
    my_genome <- BSgenome.Mmusculus.UCSC.mm10
}

library(BSgenome)
library(ensembldb)
library(stringr)
library(dplyr)

sample_id <- snakemake@wildcards$sample_id

bam.files <- list(snakemake@input$bam)
names(bam.files) <- sample_id

modelfile <- snakemake@output$modelfile

ensembldb_sqlite <- snakemake@input$ensdb
txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)

selected.genes.and.transcripts <- read.csv(snakemake@input$transcripts, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]

library(GenomicRanges)
library(alpine)

gene.lengths <- sum(width(ebt.fit))

readlength <- snakemake@params$readlength 
minsize <- 70
maxsize <- 600
gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

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


# Specify the models
models <- list(
  "all" = list(
    formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk):genelen  + GC40.80 + GC20.80 + gene",
    offset=c("fraglen","vlmm")
  )
)

# DEBUG
#trace(fitBiasModels, at=c(18,4,8), browser)
#trace(readGAlignAlpine, browser)
#trace(fitBiasModels, browser)

# Fit the models
print("Starting model fit")
fitpar <- lapply(bam.files, function(bf) {
    fitBiasModels(
        genes=ebt.fit,
        bam.file=bf,
        fragtypes=fragtypes,
        genome=my_genome,
        models=models,
        readlength=readlength,
        minsize=minsize,
        maxsize=maxsize
    )
})
print("Done model fit")

save(fitpar, file=modelfile)
