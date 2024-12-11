import pathlib
import polars as pl

fastq_dir = pathlib.Path(snakemake.input.fastq)
N_BASES = 6
MAX_NUM_SEQUENCES = 1_000_000

def zeros():
    return [0 for _ in range(N_BASES)]
def base_counts():
    return {base: zeros() for base in "ACGTN"}

forward = base_counts()
reverse = base_counts()
unpaired = base_counts()
for fastq in fastq_dir.glob("*.fastq"):
    fastq_name = fastq.name
    if "R1" in fastq_name:
        target = forward
    elif "R2" in fastq_name:
        target = reverse
    else:
        target = unpaired

    with fastq.open() as fastq:
        for _ in range(MAX_NUM_SEQUENCES):
            header = fastq.readline()
            seq = fastq.readline()
            sep = fastq.readline()
            qual = fastq.readline()

            if not header:
                break

            for i in range(N_BASES):
                target[seq[i]][i] += 1

pl.concat((
    pl.from_dicts(forward).with_columns(read=pl.lit("forward"), pos = pl.arange(N_BASES)),
    pl.from_dicts(reverse).with_columns(read=pl.lit("reverse"), pos = pl.arange(N_BASES)),
    pl.from_dicts(unpaired).with_columns(read=pl.lit("unpaired"), pos = pl.arange(N_BASES)),
)).write_csv(
    snakemake.output.hexamer,
    separator="\t"
)
