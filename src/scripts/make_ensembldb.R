library(ensembldb)

if (snakemake@wildcards$genome_name == "Mus_musculus.GRCm38.102") {
    genomeVersion <- "GRCm38.p6"
} else {
    genomeVersion <- "GRCh38.p14"
}
genomeVersion <- 
ensDbFromGtf(
     snakemake@input$gtf,
     snakemake@output$sqlite,
     genomeVersion = genomeVersion
)
