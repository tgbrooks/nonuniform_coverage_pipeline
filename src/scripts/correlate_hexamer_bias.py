import polars as pl
import scipy.ndimage
import numpy as np

#sample_id = "SRX4080514"
##sample_id = "SRX4393369"
#hexamer_file = f"data/{sample_id}/hexamer_bias.txt"
#cov_file = f"data/{sample_id}/transcript_coverage.txt"
#sequence_file = "results/selected_genomic_sequences.txt"
#strandedness = "reverse"
#outfile = f"data/{sample_id}/hexamer_correlation.txt"

sample_id = snakemake.wildcards.sample_id
hexamer_file = snakemake.input.hexamer
cov_file = snakemake.input.cov
sequence_file = snakemake.input.sequence
strandedness = snakemake.params.strandedness
outfile = snakemake.output.corr

N_PERMUTATIONS = 100

assert strandedness != "unstranded", "This script does not support unstranded data"
assert strandedness != "forward", "This script has not yet been verified on forward stranded data"

all_cov = pl.read_csv(cov_file, separator="\t")
transcript_ids = all_cov['gene'].unique()

# Load hexamer biases
raw_hexamer = pl.read_csv(hexamer_file, separator="\t") \
    .unpivot(
        on = ["A", "C", "G", "T", "N"],
        index = ["read", "pos"],
        variable_name = "base",
        value_name = "count"
    ).filter(
        pl.col('base') != 'N'
    ).with_columns(
        freq = pl.col("count") / pl.col("count").sum().over("read", "pos")
    )
hexamer = dict()
for (read,), data in raw_hexamer.group_by("read"):
    encoded = np.array([
            data.filter(base="A")['freq'],
            data.filter(base="C")['freq'],
            data.filter(base="G")['freq'],
            data.filter(base="T")['freq'],
        ],
        dtype="float32"
    )
    hexamer[read] = encoded
N_BASES = hexamer['forward'].shape[1]

def base_idx(i, complement):
    # indexes are 0 = A, 1 = C, 2 = G, 3 = T
    if complement:
        return 3 - i
    else:
        return i

# Load sequences
raw_sequences = pl.read_csv(sequence_file, separator="\t") \
    .filter(
        pl.col("transcript_id").is_in(transcript_ids)
    )

# Convert to one-hot encoded
sequences = dict()
for transcript_id, sequence in raw_sequences.iter_rows():
    x = np.frombuffer(sequence.encode("ascii"), dtype="uint8").copy()
    encoded = np.array(
        [
            x == ord("A"),
            x == ord("C"),
            x == ord("G"),
            x == ord("T"),
        ],
        dtype="float32",
    )
    sequences[transcript_id] = encoded

# Correlate to expression
def compute_hexamer_corr(hexamer):
    results = []
    for transcript_id in transcript_ids:
        cov = all_cov.filter(gene = transcript_id)['cov'].to_numpy()
        smooth_cov = scipy.ndimage.gaussian_filter1d(cov, sigma = 30, radius = 50, mode="constant", cval=0)
        small_variations = cov - smooth_cov
        diff = np.diff(small_variations)
        seq = sequences[transcript_id]

        weights = {}
        for read in ['forward', 'reverse']:
            complement = (strandedness != read)
            direction = -1 if complement else 1
            w = np.array([
                np.lib.stride_tricks.sliding_window_view(seq[i], N_BASES) * hexamer[read][base_idx(i, complement), ::direction]
                    for i in range(4)
            ]).sum(axis=0).prod(axis=1) / 0.25**N_BASES
            weights[read] = w

        if strandedness == 'reverse':
            REV_OFFSET = N_BASES
            DIFF_OFFSET = N_BASES - 1
            weight = weights['reverse'][REV_OFFSET:] - weights['forward'][:-REV_OFFSET]
            corr = np.corrcoef(diff[DIFF_OFFSET:-DIFF_OFFSET], weight)[0,1]
        else:
            raise Exception("only support reverse strandedness")

        results.append({
            "transcript_id": transcript_id,
            "corr": corr,
            })
    results = pl.DataFrame(results).with_columns(sample_id=pl.lit(sample_id))
    return results

main_results = compute_hexamer_corr(hexamer)


# Run permutations with permuted hexamer
rng = np.random.default_rng(seed=0)
def scramble_hexamer(hexamer):
    new_hexamer = dict()
    for read, hex in hexamer.items():
        new_order = rng.choice(N_BASES, size=N_BASES, replace=False)
        new_hexamer[read] = hex[:,new_order]
    return new_hexamer
scrambled = []
for i in range(N_PERMUTATIONS):
    new_hex = scramble_hexamer(hexamer)
    scrambled.append(compute_hexamer_corr(new_hex).with_columns(iter=i))
scrambled = pl.concat(scrambled)

pvalues = scrambled.join(main_results.select('transcript_id', actual_corr='corr'), on=["transcript_id"])\
    .group_by('transcript_id')\
    .agg(
        actual_corr = pl.col('actual_corr').first(),
        p = ((pl.col('corr') >= pl.col('actual_corr')).sum()+1) / (N_PERMUTATIONS+1),
        mean_perm_corr = pl.col('corr').mean(),
        std_perm_corr = pl.col('corr').std(),
        max_perm_corr = pl.col('corr').max(),
    )

pvalues.write_csv( outfile, separator="\t")
