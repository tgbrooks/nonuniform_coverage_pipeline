import pathlib
import sqlite3
import numpy as np
import polars as pl
import matplotlib as mpl
import matplotlib.pyplot as pyplot

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

sample_ids = snakemake.params.sample_ids

dupe_rate_files = snakemake.input.dupe_rate
temp = []
for (dupe_rate_file, sample_id) in zip(dupe_rate_files, sample_ids):
    dupe_rate = pl.read_csv(dupe_rate_file, has_header=False, new_columns=["chrom", "start", "end", "value"], separator="\t", dtypes=[pl.Utf8, pl.Int32, pl.Int32, pl.Float32]) \
        .with_columns(sample_id = pl.lit(sample_id))
    temp.append(dupe_rate)
dupe_rate = pl.concat(temp)

cov_files = snakemake.input.coverage
temp = []
for (cov_file, sample_id) in zip(cov_files, sample_ids):
    cov = pl.read_csv(cov_file, has_header=False, new_columns=["chrom", "start", "end", "value"], separator="\t", dtypes=[pl.Utf8, pl.Int32, pl.Int32, pl.Int32]) \
        .with_columns(sample_id = pl.lit(sample_id))
    temp.append(cov)
cov = pl.concat(temp)


annotation_db = snakemake.input.annotation
with sqlite3.connect(annotation_db) as conn:
    exons = pl.read_database("SELECT * FROM exon", conn)
    tx2exon = pl.read_database("SELECT * FROM tx2exon", conn)
    gene = pl.read_database("SELECT * FROM gene", conn)
    tx = pl.read_database("SELECT * FROM tx", conn)


# compile the list of transcripts we'll plot 
transcripts_file = snakemake.input.transcripts
## list in all the mitochondrial transcripts
mt_genes = gene.filter(pl.col("gene_name").str.starts_with("mt-"))['gene_id']
mt_tx = tx.filter(pl.col("gene_id").is_in(mt_genes))
transcripts = pl.concat((
        # Hand-picked genes
        pl.DataFrame({
            "gene_id": ["ENSMUSG00000064339","ENSMUSG00000032554", "ENSMUSG00000002985", "ENSMUSG00000072849", "ENSMUSG00000022868", "ENSMUSG00000028001", "ENSMUSG00000064337", "ENSMUSG00000064370", "ENSMUSG00000064345"],
            "transcript_id": ["ENSMUST00000082390", "ENSMUST00000112645", "ENSMUST00000174064", "ENSMUST00000085054", "ENSMUST00000023583", "ENSMUST00000166581", "ENSMUST00000082388", "ENSMUST00000082421", "ENSMUST00000082396"]
        }),
        # mitochondrial transcripts
        mt_tx.select("gene_id", pl.col("tx_id").alias("transcript_id")),
        # Our standard set of 100 single-isoform genes
        pl.read_csv(transcripts_file, separator="\t").select("gene_id", "transcript_id"),
)).filter(pl.col("gene_id").is_first_distinct())

for (gene_id, tx_id) in transcripts.iter_rows():
    chrom = gene.filter(gene_id = gene_id)['seq_name'][0]
    gene_name = gene.filter(gene_id = gene_id)['gene_name'][0]
    these_exons = (
        exons.join(tx2exon.filter(tx_id = tx_id), on="exon_id")
            .sort("exon_seq_start")
    )

    # Base-pair positions of the gene's exons in genomic coordinates
    pos = pl.concat([
        pl.DataFrame({
            "pos": np.arange(start-1, end, dtype="int32"),
            "chrom": chrom,
        }) for (start, end) in these_exons.select("exon_seq_start", "exon_seq_end").iter_rows()
    ])

    # Make the plot
    fig, axes = pyplot.subplots(figsize=(5,1.5*len(sample_ids)), nrows=len(sample_ids), layout="constrained", sharex=True, sharey=False)

    for (ax, ((sample_id,), this_cov)) in zip(axes, cov.group_by(["sample_id"], maintain_order=True)):
        this_dupe_rate = dupe_rate.filter(pl.col("sample_id") == sample_id)
        # Find the read depth/dupe_rates for each position
        # note that the cov/dupe_rate files give a range (start-end)
        # and so we have to join_asof to find the right range
        this_cov = (pos
            .join_asof(
                this_cov,
                left_on = "pos",
                right_on = "start",
                by = "chrom",
            )
            .filter(
                pl.col('pos') < pl.col('end'),
            )
        )
        this_dupe_rate = (pos
            .join_asof(
                this_dupe_rate,
                left_on = "pos",
                right_on = "start",
                by = "chrom"
            )
            .filter(
                pl.col('pos') < pl.col('end'),
            )
        )

        both = (
            this_cov.select('pos', pl.col('value').alias("cov"))
                .join(this_dupe_rate.select('pos', pl.col('value').alias('dupe_rate')), on="pos", validate="1:1")
                .with_row_index(name="loc")
                .with_columns(mask = pl.col('cov') > 100)
        )

        mask = np.array(both['mask']).astype(bool)
        masked_dupe_rate = np.array(both['dupe_rate']).astype(float)
        masked_dupe_rate[~mask] = float("nan")
        dupe_rate_lines, = ax.plot(both['loc'], masked_dupe_rate, color="b", label="dupe rate")
        ax.set_ylabel(f"{sample_id}\nPCR dupe rate")
        ax.set_xlabel("loc")
        ax_true = ax.twinx()
        coverage_lines, = ax_true.plot(both['loc'], both['cov'], color="k", label="coverage")
        ax_true.set_ylabel("coverage")
    fig.legend(handles=[dupe_rate_lines, coverage_lines], loc = "upper right")
    fig.suptitle(f"{gene_id} - {gene_name}")
    fig.savefig(outdir / f"{gene_id}.{gene_name}.png", dpi=400)
    pyplot.close()
