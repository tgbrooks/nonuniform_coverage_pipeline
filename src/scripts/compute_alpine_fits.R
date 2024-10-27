library(BSgenome)
library(ensembldb)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(splines)

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
outfile <- snakemake@output$outfile

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

# Load the BSgenome for our species
library(snakemake@params$BSgenome)
if (snakemake@params$BSgenome == "BSgenome.Hsapiens.UCSC.hg38") {
    my_genome <- BSgenome.Hsapiens.UCSC.hg38
} else {
    my_genome <- BSgenome.Mmusculus.UCSC.mm10
}


print("Making fragtypes")
make_fragtypes <- function(readlength, minsize, maxsize) {
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
  return(fragtypes)
}

# Compute the fragtypes for each model parameter value
# Since these are for different studies, they vary in some parameters that
# would normally have been held constant (read length)
get_param_str <- function(model) {
    return(paste(
        model[['model.params']][['readlength']],
        model[['model.params']][['minsize']],
        model[['model.params']][['maxsize']],
        sep=","
    ))
}
all_fragtypes <- list()
for (model in models) {
    param_str <- get_param_str(model)
    if (!(param_str %in% names(all_fragtypes))) {
      all_fragtypes[[param_str]] <- make_fragtypes(
        model[['model.params']][['readlength']],
        model[['model.params']][['minsize']],
        model[['model.params']][['maxsize']]
      )
    }
}

## Predict function - adapted from alpine
predictCoverage <- function(gene, all_fragtypes, bam.files, fitpar, genome, model.names) {
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

      param_str <- get_param_str(fitpar[[sample_id]])
      fragtypes <- all_fragtypes[[param_str]]
      fragtypes.temp <- alpine:::matchReadsToFraglist(reads, list(fragtypes))[[1]]
      ## -- fragment bias --
      fraglen.density <- fitpar[[sample_id]][["fraglen.density"]]
      stopifnot(!is.null(fraglen.density))
      fragtypes.temp$logdfraglen <- log(alpine:::matchToDensity(fragtypes.temp$fraglen,
                                                                fraglen.density))
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
    print(paste0("Plotting ", gene))
    pred.cov <- predictCoverage(gene=ebt.fit[[gene]],
                              all_fragtypes=lapply(all_fragtypes, function(fragtypes) {fragtypes[[gene]]}),
                              bam.files=bam.files,
                              fitpar=models,
                              genome=my_genome,
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
}
all_cov_table <- bind_rows(all_cov_tables)
write_delim(all_cov_table, outfile, delim="\t")
