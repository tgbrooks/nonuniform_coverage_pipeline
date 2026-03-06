import sqlite3
import pathlib
import polars as pl

AVERAGE_FRAGMENT_SIZE = 200  # just a guess
JUNCTION_RATIO_THRESHOLD = 0.05
MIN_JUNCTION_READS = 10
ensdb = snakemake.input.ensdb
quant_files = snakemake.input.quants
junction_files = snakemake.input.junctions
outfile = snakemake.output.outfile

# TESTING VALUES:
# ensdb = "data/Mus_musculus.GRCm38.102.gtf.sqlite"
# sample_ids = [
#    "SRX4080514",
#    "SRX4080520",
#    "SRX4393368",
#    "SRX4393369",
#    "SRX16386863",
#    "SRX16386864",
#    "SRX14468350",
#    "SRX14468347",
#    "SRX13396189",
#    "SRX13396186",
#    "SRX11694510",
#    "SRX11694499",
# ]
# quant_files = [f"data/{sample_id}/bam/ReadsPerGene.out.tab" for sample_id in sample_ids]
# junction_files = [f"data/{sample_id}/bam/SJ.out.tab" for sample_id in sample_ids]
# outfile = None

# Load transcript info
with sqlite3.connect(ensdb) as conn:
    gene = pl.read_database("select * from gene", conn)
    tx = pl.read_database("select * from tx", conn)
    tx2exon = pl.read_database("select * from tx2exon", conn)
    exon = pl.read_database("select * from exon", conn)

# Load gene-level quants
quants = []
for quant_file in quant_files:
    quants.append(
        pl.read_csv(
            quant_file,
            separator="\t",
            has_header=False,
            new_columns=["gene_id", "unstranded", "forward", "reverse"],
        ).with_columns(
            sample_id=pl.lit(pathlib.Path(quant_file).parent.parent.name),
        )
    )
quants = pl.concat(quants)
min_quant = quants.group_by("gene_id").agg(pl.col("unstranded").min().alias("quant"))
print(min_quant)

single_isoform = (
    tx.group_by("gene_id")
    .agg(pl.col("tx_id").count())
    .filter(pl.col("tx_id") == 1)
    .select("gene_id")
    .join(tx.select("gene_id", "tx_id"), on="gene_id")
    .join(gene, on="gene_id")
    .filter(pl.col("gene_biotype").is_in(["protein_coding", "lincRNA"]))
)


# Compute transcript lengths
# and filter
tx_lengths = (
    single_isoform.join(
        tx2exon,
        on="tx_id",
    )
    .join(
        exon,
        on="exon_id",
    )
    .group_by(["gene_id", "tx_id"])
    .agg(
        tx_length=(pl.col("exon_seq_end") - pl.col("exon_seq_start") + 1).sum(),
    )
    .filter(
        (pl.col("tx_length") > 600) & (pl.col("tx_length") < 7000),
    )
)

selected_isoforms = tx_lengths.join(single_isoform, on="tx_id")

# Check for unannotated splice junctions occuring in these transcripts
# load splice junction info from STAR's SJ.out.tab files
SJs = []
for junction_file in junction_files:
    SJs.append(
        pl.read_csv(
            junction_file,
            has_header=False,
            new_columns=[
                "chrom",
                "first_base",
                "last_base",
                "strand",
                "intron_motif",
                "annotated",
                "num_unique_reads",
                "num_multimapper_reads",
                "max_spliced_overhang",
            ],
            separator="\t",
            schema_overrides={"chrom": str},
        ).with_columns(sample_id=pl.lit(junction_file).str.split("/").list.get(1))
    )
SJ = pl.concat(SJs)

# Check for overlaps with unannotated splice junctions for our transcripts
selected_exons = (
    selected_isoforms.select("tx_id")
    .join(
        tx2exon,
        on="tx_id",
    )
    .join(
        exon,
        on="exon_id",
    )
)
problem_SJ_left = SJ.filter(
    annotated=0,
).join_where(
    selected_exons,
    pl.col("exon_seq_start") - 1 <= pl.col("first_base"),
    pl.col("exon_seq_end") + 1 >= pl.col("first_base"),
)
problem_SJ_right = SJ.filter(
    annotated=0,
).join_where(
    selected_exons,
    pl.col("exon_seq_start") - 1 <= pl.col("last_base"),
    pl.col("exon_seq_end") + 1 >= pl.col("last_base"),
)
problem_SJ = pl.concat([problem_SJ_left, problem_SJ_right]).unique()
SJ_temp = (
    problem_SJ.join(
        selected_isoforms.select("tx_id", "gene_id", "tx_length"),
        on="tx_id",
    )
    .join(
        quants.select("sample_id", "gene_id", "unstranded"), on=["gene_id", "sample_id"]
    )
    .with_columns(
        mean_cov=AVERAGE_FRAGMENT_SIZE * pl.col("unstranded") / pl.col("tx_length")
    )
    .with_columns(junction_ratio=pl.col("num_unique_reads") / pl.col("mean_cov"))
    .filter(pl.col("num_unique_reads") > MIN_JUNCTION_READS)
)
worst_SJ = SJ_temp.group_by("gene_id", "tx_id").agg(
    junction_reads=pl.col("num_unique_reads").max(),
    junction_ratio=pl.col("junction_ratio").max(),
)

high_expr = (
    selected_isoforms.join(
        min_quant,
        on="gene_id",
    )
    .join(worst_SJ.select("tx_id", "junction_ratio"), on="tx_id", how="left")
    .fill_null(0)
    .filter(pl.col("junction_ratio") < JUNCTION_RATIO_THRESHOLD)
    .sort("quant")
    .tail(100)
)
print(high_expr)
min_quant_val = high_expr.select(pl.col("quant").min())
print(f"Minimum mean expression of selected genes: {min_quant_val}")
worst_junction_ratio = high_expr.select(pl.col("junction_ratio").max())
print(
    f"Worst junction ratio (num junction reads / average coverage depth): {worst_junction_ratio}"
)

high_expr.select(
    "gene_id",
    pl.col("tx_id").alias("transcript_id"),
    "gene_name",
).write_csv(
    outfile,
    separator="\t",
)
