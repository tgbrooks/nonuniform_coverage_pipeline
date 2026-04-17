library(BSgenome)
library(ensembldb)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)
library(splines)

# Values for testing
outdir <- c("temp")
readlength <- 151
sample_id <- "SRX16386864"
tissue <- "liver"
modelfiles = paste0("data/", sample_id, "/alpine_rnafold.model.rda")
bam.files <- paste0("data/", sample_id, "/bam/Aligned.sortedByCoord.out.bam")
bai <- paste0("data/", sample_id, "/bam/Aligned.sortedByCoord.out.bam.bai")
ensembldb_sqlite <- "data/Mus_musculus.GRCm38.102.gtf.sqlite"
selected_transcripts_file <- paste0("data/", tissue, "/training_set_transcripts.txt")
fragtypes_file <- paste0("data/", tissue, "/alpine_rnafold.fragtypes.rda")
BSgenome_name <- "BSgenome.Mmusculus.UCSC.mm10"


ensembldb_sqlite <- snakemake@input$ensdb
bam.files <- snakemake@input$bam
modelfiles <- snakemake@input$models
selected_transcripts_file <- snakemake@input$transcripts
outfile <- snakemake@output$outfile
BSgenome_name <- snakemake@params$BSgenome
fragtypes_file <- snakemake@input$fragtypes
readlength <- snakemake@params$readlength

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
ebt0 <- keepSeqlevels(ebt0, standardChromosomes(ebt0), pruning.mode="coarse")
ebt.fit <- ebt0[selected.txs]

library(GenomicRanges)
library(alpine)

gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

# Load the BSgenome for our species
library(BSgenome_name, character.only=TRUE)
if (BSgenome_name == "BSgenome.Hsapiens.UCSC.hg38") {
    my_genome <- BSgenome.Hsapiens.UCSC.hg38
} else {
    my_genome <- BSgenome.Mmusculus.UCSC.mm10
}
seqlevelsStyle(my_genome) <- "NCBI"


raw_fragtypes <- readRDS(fragtypes_file)

# Fix the readlength of the fragtypes
# To avoid repeated work, we only compute fragtypes once and update its readlength to match the study
readlen_fragtypes <- sapply(
    gene.names,
    function(gene.name) {
      ft <- raw_fragtypes[[gene.name]]
      exons <- ebt.fit[[gene.name]]
      map <- alpine:::mapTxToGenome(exons)
      ft$gread1end <- alpine:::txToGenome(ft$start + readlength - 1, map)
      ft$gread2start <- alpine:::txToGenome(ft$end - readlength + 1, map)
      ft
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

get_param_str <- function(model) {
    return(paste(
        model[['model.params']][['readlength']],
        model[['model.params']][['minsize']],
        model[['model.params']][['maxsize']],
        sep=","
    ))
}

all_fragtypes <- list(readlen_fragtypes)
names(all_fragtypes) <- get_param_str(models[[sample_id]])

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
      if ((!is.null(vlmm.fivep)) && !is.null(vlmm.threep)) {
        fragtypes.temp <- alpine:::addVLMMBias(fragtypes.temp, vlmm.fivep, vlmm.threep)
      }

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
          if (models[[modeltype]]$formula == "count ~ 0 + 0") {
              models[[modeltype]]$formula <- NULL
          }
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
