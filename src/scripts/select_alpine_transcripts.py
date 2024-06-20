import sqlite3
import pathlib
import polars as pl

# Load transcript info
with sqlite3.connect(snakemake.input.ensdb) as conn:
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

# Load gene-level quants
quants = []
for quant_file in snakemake.input.quants:
    quants.append(
        pl.read_csv(
            quant_file,
            separator="\t",
            has_header = False,
            new_columns = ["gene_id", "unstranded", "forward", "reverse"],
        ).with_columns(
            id = pl.lit(pathlib.Path(quant_file).parent.parent.name),
        )
    )
quants = pl.concat(quants)
min_quant = quants.group_by("gene_id").agg(pl.col("unstranded").min().alias("quant"))
print(min_quant)

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
    .filter(
        (pl.col("tx_length") > 600) & (pl.col("tx_length") < 7000),
    )
)

selected_isoforms = tx_lengths.join(single_isoform, on="tx_id")

high_expr = (selected_isoforms
    .join(
        min_quant,
        on = "gene_id",
    )
    .sort("quant")
    .tail(100)
)
print(high_expr)
min_quant = high_expr.select(pl.col("quant").min())
print(f"Minimum mean expression of selected genes: {min_quant}")

high_expr.select(
    "gene_id",
    pl.col("tx_id").alias("transcript_id"),
    "gene_name",
).write_csv(
    snakemake.output.outfile,
    separator="\t",
)
