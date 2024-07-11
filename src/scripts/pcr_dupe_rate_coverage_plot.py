#!/usr/bin/env python

import argparse
import pysam
import numpy as np
import pandas

parser = argparse.ArgumentParser(
    "pcr_dupe_rate_coverage",
    description="Generates deduped coverage from a sorted bam file that has the 'UG' UMI group tags inserted, as in the `umit_tools group` command."
)
parser.add_argument("--bam_file", help="Sorted and indexed bam file with UG tags.")
parser.add_argument("--chrom_lengths", help="File containing two tab-separated columns with chromosome name and chromosome length")
parser.add_argument("--out_cov_file", help="File to output deduped coverage as bed format to")
parser.add_argument("--out_duped_file", help="File to output coverage as bed format to of just number of UMI groups with at least 2 copies")
parser.add_argument("--out_dupe_rate_file", help="File to output dupe rate as bed format to")

args = parser.parse_args()

START = 0
END = 1

bam_file = args.bam_file
chrom_lengths_file = args.chrom_lengths
out_cov_file = args.out_cov_file
out_duped_file = args.out_duped_file
out_dupe_rate_file = args.out_dupe_rate_file

# Clear the output files in-case, since we'll be appending to them as we go
with open(out_cov_file, "wt"):
    pass
with open(out_duped_file, "wt"):
    pass
with open(out_dupe_rate_file, "wt"):
    pass

chrom_lengths = pandas.read_csv(chrom_lengths_file, sep="\t", header=None, names=["chrom", "length"]).set_index("chrom")['length'].to_dict()

# Store one end of the read in this, mapping query id to read
partial_reads = {}

# Information for this chromosome
working_chrom = None
read_counts = np.zeros(1, np.int32)
duped_counts = np.zeros(1, np.int32)

# Map UMI groups to dupe counts
groups_encountered = set()
groups_encountered_twice = set()

def pair_blocks(first, second):
    ''' Gives the alignment blocks of a read pair, removing overlap '''
    blocks1 = first.blocks
    blocks2 = second.blocks
    assert blocks1[0][START] <= blocks2[0][START]

    if blocks1[-1][END] > blocks2[0][START]:
        # Overlap! We'll trim blocks1 to not overlap
        new_end = blocks2[0][START]
        blocks1 = [(s,min(e, new_end)) for (s,e) in blocks1
                        if s < new_end]
    yield from blocks1
    yield from blocks2

def rle(x):
    ''' Run-length encode - groups runs of the same values together '''
    breakpoints, = np.where(np.diff(x, prepend=0, append=0))
    values = x[breakpoints]
    non_zero = np.where(values)
    starts = breakpoints[:-1]
    ends = breakpoints[1:]
    return starts[non_zero], ends[non_zero], values[non_zero]

def to_bed(x, chrom):
    (s,e,v) = rle(x)
    return pandas.DataFrame({
        "chrom": chrom,
        "start": s,
        "end": e,
        "value": v,
    })

def finish_chromosome():
    print(f"Finished with {working_chrom}")
    print(f"Processed {read_counts.sum()} deduped")
    print(f"Processed {duped_counts.sum()} dupes")

    dupe_rate = np.divide(duped_counts, read_counts, out = np.zeros(len(read_counts)), where = read_counts > 0)
    with open(out_cov_file, "at") as out:
        to_bed(read_counts, working_chrom).to_csv(out, index=False, header=False, sep="\t")
    with open(out_duped_file, "at") as out:
        to_bed(duped_counts, working_chrom).to_csv(out, index=False, header=False, sep="\t")
    with open(out_dupe_rate_file, "at") as out:
        to_bed(dupe_rate, working_chrom).to_csv(out, index=False, header=False, sep="\t")

print(f"Processing {bam_file}")
i = 0
with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam:
        i += 1
        if i % 100_000 == 0:
            print(i)

        if read.is_secondary:
            continue

        try:
            other = partial_reads.pop(read.query_name)
        except KeyError:
            # We haven't processed the pair of the paired ends yet
            partial_reads[read.query_name] = read 
            continue

        # At this point, we have a read pair
        try:
            umi_group = read.get_tag("UG")
        except KeyError:
            # Each read id has exactly one entry with a UG read
            # So now try the other half of the pair
            try:
                umi_group = other.get_tag("UG")
            except KeyError:
                raise Exception(f"The read {read.query_name} had no UMI group tag")

        if umi_group not in groups_encountered:
            # First time encountering this group
            groups_encountered.add(umi_group)
            for (s,e) in pair_blocks(other, read):
                read_counts[s:e] += 1
        elif umi_group not in groups_encountered_twice:
            # This group contains a duplicate
            groups_encountered_twice.add(umi_group)
            for (s,e) in pair_blocks(other, read):
                duped_counts[s:e] += 1
        else:
            # Do nothing - we have already counted this group *and* marked it as containing a duplicate
            pass

        if working_chrom != read.reference_name:
            if working_chrom != None:
                finish_chromosome()

            # initialize the new chromosome
            working_chrom = read.reference_name
            chrom_length = chrom_lengths[working_chrom]
            read_counts = np.zeros(chrom_length+1, np.int32)
            duped_counts = np.zeros(chrom_length+1, np.int32)
            groups_encountered = set()

finish_chromosome()
print(f"{len(partial_reads)} reads had no matching pair")
