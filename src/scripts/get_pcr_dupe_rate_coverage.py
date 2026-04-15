import sqlite3
import numpy as np
import polars as pl

annotation_db = snakemake.input.annotation
transcripts_file = snakemake.input.transcripts
sample_ids = snakemake.params.sample_ids


with sqlite3.connect(annotation_db) as conn:
    exons = pl.read_database("SELECT * FROM exon", conn)
    tx2exon = pl.read_database("SELECT * FROM tx2exon", conn)
    gene = pl.read_database("SELECT * FROM gene", conn)
    tx = pl.read_database("SELECT * FROM tx", conn)


# compile the list of transcripts we'll plot
## list in all the mitochondrial transcripts
mt_genes = gene.filter(pl.col("gene_name").str.starts_with("mt-"))["gene_id"]
mt_tx = tx.filter(pl.col("gene_id").is_in(mt_genes))
transcripts = pl.concat(
    (
        ## Hand-picked genes
        # pl.DataFrame(
        #    {
        #        "gene_id": [
        #            "ENSMUSG00000064339",
        #            "ENSMUSG00000032554",
        #            "ENSMUSG00000002985",
        #            "ENSMUSG00000072849",
        #            "ENSMUSG00000022868",
        #            "ENSMUSG00000028001",
        #            "ENSMUSG00000064337",
        #            "ENSMUSG00000064370",
        #            "ENSMUSG00000064345",
        #        ],
        #        "transcript_id": [
        #            "ENSMUST00000082390",
        #            "ENSMUST00000112645",
        #            "ENSMUST00000174064",
        #            "ENSMUST00000085054",
        #            "ENSMUST00000023583",
        #            "ENSMUST00000166581",
        #            "ENSMUST00000082388",
        #            "ENSMUST00000082421",
        #            "ENSMUST00000082396",
        #        ],
        #    }
        # ),
        ## mitochondrial transcripts
        # mt_tx.select("gene_id", pl.col("tx_id").alias("transcript_id")),
        # Our standard set of 100 single-isoform genes
        pl.read_csv(transcripts_file, separator="\t").select(
            "gene_id", "transcript_id"
        ),
    )
).filter(pl.col("gene_id").is_first_distinct())


def get_value(pos, cov):
    # Maps positions to an RLE-encoded coverage ala bed files
    # assumes pos contains one gene on one chromosome
    gene_end = pos["pos"].max()
    gene_start = pos["pos"].min()
    chrom = pos["chrom"][0]
    small_cov = (
        cov.filter(  # this turns out to much faster than letter Polars do the full join
            pl.col("start") <= gene_end,
            pl.col("end") >= gene_start,
            chrom=chrom,
        )
    ).sort("start")
    return (
        pos.join_asof(
            small_cov,
            left_on="pos",
            right_on="start",
        )
        .with_columns(pl.col("value").fill_null(0))
        .with_columns(
            value=pl.when(pl.col("end") > pl.col("pos"))
            .then(pl.col("value"))
            .otherwise(pl.lit(0))
        )
        .select("pos", "chrom", "value")
    )


def get_cov_and_dupe(gene_id, cov, dupe_rate, starts, duped_starts):
    tx_id = str(transcripts.filter(gene_id=gene_id)["transcript_id"][0])
    chrom = gene.filter(gene_id=gene_id)["seq_name"][0]
    these_exons = exons.join(tx2exon.filter(tx_id=tx_id), on="exon_id").sort(
        "exon_seq_start"
    )

    # Base-pair positions of the gene's exons in genomic coordinates
    pos = pl.concat(
        [
            pl.DataFrame(
                {
                    "pos": np.arange(start - 1, end, dtype="int32"),
                    "chrom": chrom,
                }
            )
            for (start, end) in these_exons.select(
                "exon_seq_start", "exon_seq_end"
            ).iter_rows()
        ]
    ).set_sorted("pos")

    # Find the read depth/dupe_rates for each position
    # note that the cov/dupe_rate files give a range (start-end)
    # and so we have to join_asof to find the right range
    this_cov = get_value(pos, cov)
    this_dupe_rate = get_value(pos, dupe_rate)
    this_starts = get_value(pos, starts)
    this_duped_starts = get_value(pos, duped_starts)

    both = (
        this_cov.select("pos", pl.col("value").alias("cov"))
        .join(
            this_dupe_rate.select("pos", pl.col("value").alias("dupe_rate")),
            on="pos",
            validate="1:1",
        )
        .join(
            this_starts.select("pos", pl.col("value").alias("starts")),
            on="pos",
            validate="1:1",
        )
        .join(
            this_duped_starts.select("pos", pl.col("value").alias("duped_starts")),
            on="pos",
            validate="1:1",
        )
        .with_row_count(name="loc")
        .with_columns(
            transcript_id=pl.lit(tx_id),
        )
    )
    return both


def read_bed(file):
    try:
        return pl.read_csv(
            file,
            has_header=False,
            new_columns=["chrom", "start", "end", "value"],
            separator="\t",
            dtypes=[pl.Utf8, pl.Int32, pl.Int32, pl.Float32],
        )
    except pl.NoDataError:
        # Polars errors if a csv is empty
        # despite it already having the full schema
        return pl.DataFrame(
            data=[],
            schema={
                "chrom": pl.Utf8,
                "start": pl.Int32,
                "end": pl.Int32,
                "value": pl.Float32,
            },
        )


# Load the data and transform it into a single dataframe
temp = []
for sample_id in sample_ids:
    print(sample_id)
    coverage_file = (
        f"data/{sample_id}/pcr_dupe_rate_coverage/{sample_id}.deduped_coverage.bed"
    )
    dupe_rate_file = f"data/{sample_id}/pcr_dupe_rate_coverage/{sample_id}.pcr_dupe_rate_coverage.bed"
    starts_file = (
        f"data/{sample_id}/pcr_dupe_rate_coverage/{sample_id}.deduped_starts.bed"
    )
    duped_starts_file = (
        f"data/{sample_id}/pcr_dupe_rate_coverage/{sample_id}.duped_starts.bed"
    )

    all_dupe_rate = read_bed(dupe_rate_file)
    all_coverage = read_bed(coverage_file)
    all_starts = read_bed(starts_file)
    all_duped_starts = read_bed(duped_starts_file)

    for gene_id in transcripts["gene_id"]:
        temp.append(
            get_cov_and_dupe(
                gene_id,
                all_coverage,
                all_dupe_rate,
                all_starts,
                all_duped_starts,
            ).with_columns(
                sample_id=pl.lit(sample_id),
                gene_id=pl.lit(gene_id),
            )
        )
data = pl.concat(temp)
summed_data = (
    data.group_by("transcript_id", "pos", "loc")
    .agg(
        cov=pl.col("cov").sum(),
        # combine dupe rates across the different samples, weighted sum according to the depth of coverage
        dupe_rate=(
            (pl.col("cov") * pl.col("dupe_rate")).sum() / pl.col("cov").sum()
        ).fill_nan(pl.lit(None)),
        starts=pl.col("starts").sum(),
        duped_starts=pl.col("duped_starts").sum(),
    )
    .sort(["transcript_id", "loc"])
)

data.write_csv(snakemake.output.full_cov, separator="\t")
summed_data.write_csv(snakemake.output.summed_cov, separator="\t")
