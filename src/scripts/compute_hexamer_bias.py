import pathlib
import polars as pl
import pysam

bamfile = snakemake.input.bam
#bamfile = "data/SRX4393368/bam/Aligned.sortedByCoord.out.bam"
#bamfile = "data/SRX16386864/bam/Aligned.sortedByCoord.out.bam"

N_BASES = 12
MAX_NUM_SEQUENCES = 500_000_000

def zeros():
    return [0 for _ in range(N_BASES)]
def base_counts():
    return {base: zeros() for base in "ACGTN"}

forward = base_counts()
reverse = base_counts()

match_ops = set([pysam.CMATCH, pysam.CINS, pysam.CEQUAL, pysam.CDIFF])
consume_ops = set([pysam.CMATCH, pysam.CINS, pysam.CEQUAL, pysam.CDIFF, pysam.CSOFT_CLIP, pysam.CHARD_CLIP])
def mapped_bases(read):
    # Indexes of bases in the query that 'count' for hexamer priming
    # in particular, don't count soft/hard clipped bases
    cigartuple = read.cigartuples
    if read.is_reverse:
        cigartuple = cigartuple[::-1]
    start = 0
    for (op, length) in cigartuple:
        if op in match_ops:
            yield from range(start, start+length)
        if op in consume_ops:
            start += length

samfile = pysam.AlignmentFile(bamfile, 'rb')
for i, read in enumerate(samfile.fetch()):
    if i > MAX_NUM_SEQUENCES:
        break

    if read.is_secondary:
        continue

    if read.is_read1:
        target = forward
    else:
        target = reverse

    seq = read.get_forward_sequence()
    for j, query_idx in enumerate(mapped_bases(read)):
        if j >= N_BASES:
            break
        letter = seq[query_idx]
        target[letter][j] += 1

pl.concat((
    pl.from_dicts(forward).with_columns(read=pl.lit("forward"), pos = pl.arange(N_BASES)),
    pl.from_dicts(reverse).with_columns(read=pl.lit("reverse"), pos = pl.arange(N_BASES)),
)).write_csv(
    snakemake.output.hexamer,
    separator="\t"
)
