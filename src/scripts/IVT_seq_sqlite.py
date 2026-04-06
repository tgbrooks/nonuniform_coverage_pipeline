import sqlite3
import subprocess
import io
import pathlib
import requests
import gzip
import re
import polars as pl

# Download the original annotations from the IVT-seq paper
annot_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1219398&format=file&file=GSM1219398%5Fivt%2Donly%5Frep1%5Ffeature%5Fquant%2Etxt%2Egz"
response = requests.get(annot_URL)
response.raise_for_status()
annot_raw = gzip.decompress(response.content).decode()

# We also need to convert from hg19 to GRCh38 since the IVT_seq annotations are in hg19
chain_file_url = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)
chain_file_path = pathlib.Path("temp.hg19ToHg38.over.chain.gz")
if not chain_file_path.exists():
    chain_file_response = requests.get(chain_file_url)
    chain_file_response.raise_for_status()
    with open(chain_file_path, "wb") as f:
        f.write(chain_file_response.content)

ensdb = "data/GRCh38.ensemblv109.gtf.sqlite"
with sqlite3.connect(ensdb) as conn:
    gene = pl.read_database("select * from gene", conn)
    # tx = pl.read_database("select * from tx", conn)
    # tx2exon = pl.read_database("select * from tx2exon", conn)
    # exon = pl.read_database("select * from exon", conn)


def chrom_name(chrom):
    if chrom == "chrM":
        return "MT"
    else:
        return chrom.removeprefix("chr")


def parse_exons(annot_raw):
    rows = []
    current_gene = None
    current_strand = None

    for line in annot_raw.splitlines():
        line = line.rstrip()

        # Match gene/transcript header like "BC008991(GenBank)	-"
        m = re.match(r"^(\S+)\((\S+)\)\s+([-+])", line)
        if m:
            current_gene = m.group(1)
            current_strand = {"-": -1, "+": 1}[m.group(3)]
            continue

        # Match exon lines
        m = re.match(r"^\s+exon\s+(\d+)\s+(\w+):(\d+)-(\d+)\s", line)
        if m and current_gene:
            rows.append(
                {
                    "transcript_id": current_gene,
                    "exon_number": int(m.group(1)),
                    "chrom": m.group(2),
                    "start": int(m.group(3)),
                    "end": int(m.group(4)),
                    "strand": current_strand,
                }
            )

    return pl.DataFrame(rows)


hg19_annot = parse_exons(annot_raw).with_columns(
    exon_id=pl.col("transcript_id") + "_exon" + pl.col("exon_number").cast(str),
)

# Carry-over to GRCh38
bed_file = pathlib.Path("temp.annot.bed")
hg19_annot.select("chrom", "start", "end").write_csv(
    bed_file, separator="\t", include_header=False
)
res = subprocess.run(
    f"CrossMap bed {chain_file_path} {bed_file}",
    shell=True,
    check=True,
    stdout=subprocess.PIPE,
)
chained_annot = pl.read_csv(
    io.StringIO(res.stdout.decode()),
    separator="\t",
    has_header=False,
    new_columns=[
        "orig_chrom",
        "orig_start",
        "orig_end",
        "sep",
        "new_chrom",
        "new_start",
        "new_end",
    ],
)

GRCh38_annot = hg19_annot.join(
    chained_annot,
    left_on=["chrom", "start", "end"],
    right_on=["orig_chrom", "orig_start", "orig_end"],
    how="left",
).select(
    "transcript_id",
    pl.col("chrom").map_elements(chrom_name, return_dtype=str),
    start="new_start",
    end="new_end",
    strand="strand",
    exon_number="exon_number",
    exon_id="exon_id",
)

# Annotate the ENSEMBL genes
matches = GRCh38_annot.join_where(
    gene,
    pl.col("chrom") == pl.col("seq_name"),
    pl.col("strand") == pl.col("seq_strand"),
    pl.col("end") >= pl.col("gene_seq_start"),
    pl.col("start") <= pl.col("gene_seq_end"),
)
gene_matches = (
    matches.filter(gene_biotype="protein_coding")
    .select("transcript_id", "gene_id", "gene_name")  # , "gene_biotype")
    .sort("gene_id")
    .group_by("transcript_id")
    .agg(
        pl.col("gene_id").first(),
        pl.col("gene_name").first(),
    )
)

annot = GRCh38_annot.join(
    gene_matches,
    "transcript_id",
    how="left",
)


ivt_gene = annot.select(
    gene_id="gene_id",
    gene_name="gene_name",
    seq_strand="strand",
    seq_name="chrom",
    gene_biotype=pl.lit("protein_coding"),
).unique()
ivt_exon = annot.select("exon_id", exon_seq_start="start", exon_seq_end="end").unique()
ivt_tx = annot.select(
    tx_id="transcript_id",
    gene_id="gene_id",
    tx_biotype=pl.lit("protein_coding"),
).unique()
ivt_tx2exon = annot.select(
    tx_id="transcript_id", exon_id="exon_id", exon_idx="exon_number"
).unique()


# We use this function to avoid the pyarrow/sqlalchemy/pandas requirements
def write_to_sqlite(df, table_name, db_path, schema: dict | None = None):
    conn = sqlite3.connect(db_path)
    conn.execute(f"DROP TABLE IF EXISTS {table_name}")
    types = schema or {c: "TEXT" for c in df.columns}
    cols = ", ".join(f"{c} {types[c]}" for c in df.columns)
    conn.execute(f"CREATE TABLE {table_name} ({cols})")
    conn.executemany(
        f"INSERT INTO {table_name} VALUES ({','.join('?' * len(df.columns))})",
        df.rows(),
    )
    conn.commit()
    conn.close()


out_db = "data/IVT_seq.gtf.sqlite"
write_to_sqlite(
    ivt_gene,
    "gene",
    out_db,
    schema={
        "gene_id": "TEXT",
        "seq_strand": "INTEGER",
        "seq_name": "TEXT",
        "gene_name": "TEXT",
        "gene_biotype": "TEXT",
    },
)
write_to_sqlite(
    ivt_exon,
    "exon",
    out_db,
    schema=dict(exon_id="TEXT", exon_seq_start="INTEGER", exon_seq_end="INTEGER"),
)
write_to_sqlite(ivt_tx, "tx", out_db)
write_to_sqlite(
    ivt_tx2exon,
    "tx2exon",
    out_db,
    schema=dict(tx_id="TEXT", exon_id="TEXT", exon_idx="INTEGER"),
)


# Clean up temp files
chain_file_path.unlink()
bed_file.unlink()
