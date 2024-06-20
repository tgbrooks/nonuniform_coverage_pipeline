#options(error = function() traceback(10))
options(error = function() recover())

outdir <- snakemake@output$outdir
if (!dir.exists(outdir)) {dir.create(outdir)}

library(alpine)

#modelfiles <- c("data/SRX4080514/alpine.model.rda", "data/SRX4393368/alpine.model.rda", "data/SRX16386863/alpine.model.rda")
modelfiles <- snakemake@input$models
print(modelfiles)
models <- lapply(
    modelfiles,
    function (modelfile) {
      load(modelfile)
      return(fitpar)
    }
)
models <- unlist(models, recursive=FALSE)
print("Loaded:")
print(names(models))

library(BSgenome)
library(ensembldb)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(splines)

#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
ensembldb_sqlite <- snakemake@input$ensdb
txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]

selected.genes.and.transcripts <- read.csv(snakemake@input$transcripts, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]
print(head(ebt.fit))

library(GenomicRanges)
library(alpine)

### Plot functions - adapted from alpine
plotRelPos <- function(fitpar, model, col, lty, ylim, genelen) {
  # Plot rel pos at a given gene length since we include a genelen interaction term

  # just a single sample?
  if ("models" %in% names(fitpar)) {
    fitpar <- list(fitpar)
  }

  knots <- fitpar[[1]][["model.params"]]$relpos.knots
  bk <- fitpar[[1]][["model.params"]]$relpos.bk

  all_fitpars <- fitpar
  plot(0,0,type="n",xlim=c(0,1),
      ylab="log fragment rate", xlab="5' -- position in transcript -- 3'",
      main=paste0("relative position bias (", genelen / 1000, "kb)"))
  for (this_fitpar in names(all_fitpars)) {
    fitpar <- all_fitpars[this_fitpar]
    n <- length(knots)
    coef.nms <- names(fitpar[[1]][["coefs"]][[model]])
    coef.idx <- c(grep("\\(Intercept\\)",coef.nms), grep("ns\\(relpos", coef.nms))
    coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[model]][coef.idx])
    z <- seq(from=0,to=1,length=101)
    genelen_ <- rep(genelen, length(z))
    # include the genelength in the model
    x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk) + ns(z, knots=knots, Boundary.knots=bk):genelen_)
    logpred <- x %*% coefmat
    logpred <- scale(logpred, scale=FALSE)
    if (missing(ylim)) {
      ylim <- c(min(logpred),max(logpred))
    }
    if (missing(col)) {
      col <- rep("black", ncol(logpred))
    }
    if (missing(lty)) {
      lty <- rep(1, ncol(logpred))
    }
    for (i in 1:ncol(logpred)) {
      lines(z, logpred[,i], col=col[names(all_fitpars) == this_fitpar], lwd=2, lty=lty[i])
    }
  }
}

library(RColorBrewer)
ids <- names(models)
ids_coded <- as.integer(factor(ids))
color_by_id <- brewer.pal(length(ids), "Set1")

# make plots
print("Starting to generate frag length image")
palette(color_by_id)
png(paste(outdir, "frag_len.by_id.png", sep="/"),  width=6, height=6, res=400, units="in")
plotFragLen(models, col=ids_coded)
legend("topleft", legend=unique(ids),
              col=color_by_id[unique(ids_coded)],
              lty=1, lwd=4)
dev.off()

print("Starting to generate GC image")
png(paste(outdir, "GC.by_id.png", sep="/"),  width=6, height=6, res=400, units="in")
plotGC(models, model="all", col=ids_coded)
legend("bottomright", legend=unique(ids),
              col=color_by_id[unique(ids_coded)],
              lty=1, lwd=4)
dev.off()

png(paste(outdir, "RelPos.by_id.genelen=1000.png", sep="/"),  width=6, height=6, res=400, units="in")
plotRelPos(models, model="all", col=ids_coded, genelen=1000, ylim=c(-1,1))
legend("bottomright", legend=unique(ids),
              col=color_by_id[unique(ids_coded)],
              lty=1, lwd=4)
dev.off()

png(paste(outdir, "RelPos.by_id.genelen=5000.png", sep="/"),  width=6, height=6, res=400, units="in")
plotRelPos(models, model="all", col=ids_coded, genelen=5000, ylim=c(-1,1))
legend("bottomright", legend=unique(ids),
              col=color_by_id[unique(ids_coded)],
              lty=1, lwd=4)
dev.off()
