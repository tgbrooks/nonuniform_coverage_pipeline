## COMPUTE SINGLE-ISOFORM GENES IN HUMANS
## NOT CONSIDERING WHETHER HIGH-EXPRESSED

import sqlite3
import pathlib
import polars as pl


genome_name = "GRCh38.ensemblv109"
ensdb = f"data/{genome_name}.gtf.sqlite"

# Load transcript info
with sqlite3.connect(ensdb) as conn:
    gene = pl.read_database(
        "select * from gene",
        conn
    )
    tx = pl.read_database(
        "select * from tx",
        conn
    )
    tx2exon = pl.read_database(
        "select * from tx2exon",
        conn
    )
    exon = pl.read_database(
        "select * from exon",
        conn
    )

single_isoform = (tx
    .group_by("gene_id")
    .agg(pl.col("tx_id").count())
    .filter(pl.col('tx_id') == 1)
    .select('gene_id')
    .join(tx.select('gene_id', 'tx_id'), on="gene_id")
    .join(gene, on="gene_id")
    .filter(pl.col("gene_biotype").is_in(["protein_coding", "lincRNA"]))
)


# Compute transcript lengths
# and filter
tx_lengths = (single_isoform
    .join(
        tx2exon,
        on = "tx_id",
    )
    .join(
        exon,
        on = "exon_id",
    )
    .group_by(["gene_id", "tx_id"])
    .agg(
        tx_length = (pl.col("exon_seq_end") - pl.col("exon_seq_start") + 1).sum(),
    )
)

tx_lengths.write_csv("results/scratch/human_single_isoform_genes.txt", separator="\t")
