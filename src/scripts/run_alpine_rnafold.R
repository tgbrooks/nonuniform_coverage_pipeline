#options(error = function() traceback(10))

library(BSgenome)
library(ensembldb)
library(stringr)
library(r2r)

#sample_id <- "SRX16386863"
#readlength <- 151
#modelfile <- "results/alpine_rnafold/SRX16386863/model.rda"
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#bam.files <- list("data/SRX16386863/bam/Aligned.sortedByCoord.out.bam")
#transcript_file <- "data/liver/high_expressed_single_isoform_genes.txt"

RNA_FOLD_RESOLUTION <- 25

sample_id <- snakemake@wildcards$sample_id
readlength <- snakemake@params$readlength 
modelfile <- snakemake@output$modelfile
ensembldb_sqlite <- snakemake@input$ensdb
bam.files <- list(snakemake@input$bam)
transcript_file <- snakemake@input$transcripts

names(bam.files) <- sample_id

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

# Load Mmusculus object
# For our prepared BSgenome object
# See https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
# for instructions on how it was prepared
library(BSgenome.Mmusculus.ENSEMBL.GRCm38)

fragtypes <- readRDS(snakemake@input$fragtypes)

# Specify the models
models <- list(
  "all" = list(
    formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk):genelen  + GC40.80 + GC20.80 + MFE_per_nt + I(MFE_per_nt^2) + gene",
    offset=c("fraglen","vlmm")
  )
)

# Fit the models
print("Starting model fit")
fitpar <- lapply(bam.files, function(bf) {
    fitBiasModels(
        genes=ebt.fit,
        bam.file=bf,
        fragtypes=fragtypes,
        genome=Mmusculus,
        models=models,
        readlength=readlength,
        minsize=minsize,
        maxsize=maxsize
    )
})
print("Done model fit")

save(fitpar, file=modelfile)
