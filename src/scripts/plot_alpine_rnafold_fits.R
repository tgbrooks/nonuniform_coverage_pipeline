library(BSgenome)
library(ensembldb)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(splines)

options(error = quote({dump.frames(); save.image(file = "last.dump.rda"); quit(save="no", status=1)}))

# Values for testing
#outdir <- c("results/alpine_fit_plots")
#modelfiles <- c("data/SRX4080514/alpine.model.rda", "data/SRX4393368/alpine.model.rda", "data/SRX16386863/alpine.model.rda")
#bam.files <- c("data/SRX4080514/bam/Aligned.sortedByCoord.out.bam", "data/SRX4393368/bam/Aligned.sortedByCoord.out.bam", "data/SRX16386863/bam/Aligned.sortedByCoord.out.bam")
#ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
#selected_transcripts_file <- "data/liver/high_expressed_single_isoform_genes.txt"

ensembldb_sqlite <- snakemake@input$ensdb
bam.files <- snakemake@input$bam
modelfiles <- snakemake@input$models
selected_transcripts_file <- snakemake@input$transcripts
fragtypes_file <- snakemake@input$fragtypes
outdir <- snakemake@output$outdir

dir.create(outdir)

# Load the alpine models
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
sample_ids <- names(models)


names(bam.files) <- sample_ids

meta_data <- tibble(
    sample = sample_ids,
)

txdb <- EnsDb(ensembldb_sqlite)
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]

selected.genes.and.transcripts <- read.csv(selected_transcripts_file, sep="\t")
selected.txs = selected.genes.and.transcripts$transcript_id

ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[selected.txs]

library(GenomicRanges)
library(alpine)

gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

# Load Mmusculus object
# For our prepared BSgenome object
# See https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
# for instructions on how it was prepared
library(BSgenome.Mmusculus.ENSEMBL.GRCm38)

fragtypes <- readRDS(fragtypes_file)
adaptFragtypes <- function(fragtypes, readlength, gene) {
  # Fix the readlength of the fragtypes
  # To avoid repeated work, we only compute fragtypes once and update its readlength to match the study
  map <- alpine:::mapTxToGenome(gene)
  fragtypes$gread1end <- alpine:::txToGenome(fragtypes$start + readlength - 1, map)
  fragtypes$gread2start <- alpine:::txToGenome(fragtypes$end - readlength + 1, map)
  return(fragtypes)
}

## Predict function - adapted from alpine
predictCoverage <- function(gene, gene_fragtypes, bam.files, fitpar, genome, model.names) {
    stopifnot(is(gene, "GRanges"))
    stopifnot(!is.null(fitpar))
    stopifnot(all(names(bam.files) %in% names(fitpar)))
    if (is.null(names(bam.files))) {
      names(bam.files) <- seq_along(bam.files)
    }

    # take model names and fitpar models and make the
    # models suitable for bias calculation
    models <- alpine:::namesToModels(model.names, fitpar)

    res <- list()
    for (sample_id in names(bam.files)) {
      # add counts
      bam.file <- bam.files[sample_id]
      generange <- range(gene)
      strand(generange) <- "*" # not necessary
      suppressWarnings({
        ga <- alpine:::readGAlignAlpine(bam.file, generange)
      })
      if (length(ga) == 0) {
        res[[sample_id]] <- as.list(rep(NA,length(models)))
        names(res[[sample_id]]) <- names(models)
        next
      }
      ga <- GenomeInfoDb::keepSeqlevels(ga, as.character(seqnames(gene)[1]))
      fco <- GenomicAlignments::findCompatibleOverlaps(ga, GRangesList(gene))
      # message("-- ",round(length(fco)/length(ga),2)," compatible overlaps")
      reads <- alpine:::gaToReadsOnTx(ga, GRangesList(gene), fco)

      # save fragment coverage for later
      l <- sum(width(gene))
      frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])


      fragtypes.temp <- adaptFragtypes(gene_fragtypes, fitpar[[sample_id]]$model.params$readlength, gene)
      fragtypes.temp <- alpine:::matchReadsToFraglist(reads, list(fragtypes.temp))[[1]]
      ## -- fragment bias --
      fraglen.density <- fitpar[[sample_id]][["fraglen.density"]]
      stopifnot(!is.null(fraglen.density))
      fragtypes.temp$logdfraglen <- log(alpine:::matchToDensity(fragtypes.temp$fraglen, fraglen.density))
      ## -- random hexamer priming bias with VLMM --
      vlmm.fivep <- fitpar[[sample_id]][["vlmm.fivep"]]
      vlmm.threep <- fitpar[[sample_id]][["vlmm.threep"]]
      stopifnot(!is.null(vlmm.fivep))
      stopifnot(!is.null(vlmm.threep))
      fragtypes.temp <- alpine:::addVLMMBias(fragtypes.temp, vlmm.fivep, vlmm.threep)

      # -- fit models --
      res[[sample_id]] <- list()

      # remove first and last bp for predicting coverage along transcript
      not.first.or.last.bp <- !(fragtypes.temp$start == 1 | fragtypes.temp$end == l)
      fragtypes.temp <- fragtypes.temp[not.first.or.last.bp,]

      ir <- IRanges(fragtypes.temp$start, fragtypes.temp$end)
      res[[sample_id]]$l <- l
      res[[sample_id]]$frag.cov <- frag.cov
      res[[sample_id]]$pred.cov <- list()
      for (modeltype in names(models)) {
        # message("predicting model type: ",modeltype)
        log.lambda <- alpine:::getLogLambda(fragtypes.temp, models, modeltype, fitpar, sample_id)
        pred0 <- exp(log.lambda)
        pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
        res[[sample_id]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
      }
    }
    res
}

# Plot the estimates of read depths
model.names = c("all")

all_cov_tables <- list()
for( gene in names(ebt.fit) ) {
  print(paste0("Fitting ", gene))
  pred.cov <- predictCoverage(gene=ebt.fit[[gene]],
                              gene_fragtypes=fragtypes[[gene]],
                              bam.files=bam.files,
                              fitpar=models,
                              genome=Mmusculus,
                              model.names=model.names)

  # Assemble the actual+predicted data into a long-form tibble
  # for use with ggplot
  cov_table <- tibble()
  for( sample in names(pred.cov)) {
    predicted <- as.vector(pred.cov[[sample]]$pred.cov$all)
    actual <- as.vector(pred.cov[[sample]]$frag.cov)
    max_length <- max(length(predicted), length(actual), 100) # Fills in zeros to 100 bp if we had no reads
    this_cov_table <- tibble(
          sample = sample,
          gene = gene,
          pos = 1:max_length,
          predicted = c(predicted, rep(0, max_length-length(predicted))),
          actual = c(actual, rep(0, max_length-length(actual))),
    )
    cov_table <- bind_rows(cov_table, this_cov_table)
  }
  cov_table <- cov_table %>%
    left_join(
              meta_data,
              by = "sample",
    )
  all_cov_tables[[length(all_cov_tables)+1]] <- cov_table

  print(paste("Plotting", gene))
  gene_cov <- cov_table[cov_table$gene == gene,]
  ggplot(
         data = gene_cov,
         aes(x=pos, y=actual)
    ) +
    facet_wrap(
                ~sample,
                scales = "free_y",
                ncol = 1,
                ) +
    geom_path(
              color = "black",
              ) +
    geom_path(
              aes(y = predicted),
              color = "red"
              ) +
    labs(
          x = "Position",
          y = "Read Depth",
          title = "Alpine model fits"
          ) +
    theme(
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    )
  ggsave(
          paste(outdir, paste("trace.", gene, ".png", sep=""), sep="/"),
          width = 7,
          height = 7,
  )
}

write_tsv(bind_rows(all_cov_tables), paste(outdir, "all_cov.txt"))
