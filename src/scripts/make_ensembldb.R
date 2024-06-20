library(ensembldb)

ensDbFromGtf(
     snakemake@input$gtf,
     snakemake@output$sqlite,
)
