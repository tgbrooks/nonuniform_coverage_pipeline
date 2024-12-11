import polars as pl
import scipy.ndimage
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pyplot

sample_id = "SRX4080514"
hexamer_file = f"data/{sample_id}/hexamer_bias.txt"
cov_file = f"data/{sample_id}/transcript_coverage.txt"
sequence_file = "results/selected_genomic_sequences.txt"
strandedness = "reverse"

assert strandedness != "unstranded", "This script does not support unstranded data"
assert strandedness != "forward", "This script has not yet been verified on forward stranded data"

all_cov = pl.read_csv(cov_file, separator="\t")
transcript_ids = all_cov['gene'].unique()

transcript_id = "ENSMUST00000082402" # the mean gene for correlation in this sample

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

# MAKE THE FIGURE
fig, (ax, ax1, ax2) = pyplot.subplots(figsize=(6,6), sharex = True, nrows=3)
pos = np.arange(len(cov)) / 1_000
dpos1 = pos[:-1]
dpos2 = pos[DIFF_OFFSET:-DIFF_OFFSET-1]
ax.plot(
    pos,
    cov,
    color = "r",
    label = "coverage"
)
ax.plot(pos, smooth_cov, color = "k", label = "smoothed")
ax1.plot(dpos1, diff, color="grey", label="diff")
x1, x2 = 0.2,0.5 # inset bounds
inset_values = diff[(dpos1 >= x1) & (dpos1 <= x2)]
y1, y2 = inset_values.min(), inset_values.max()
axin1 = ax1.inset_axes([0.2, 0.05, 0.6, 0.25], xlim=(x1,x2), ylim=(y1,y2), xticks=[], yticks=[])
ax1.indicate_inset_zoom(axin1, edgecolor="black")
axin1.plot(dpos1, diff, color="grey")
ax2.plot(dpos2, weight)
inset_values = weight[(dpos2 >= x1) & (dpos2 <= x2)]
y1, y2 = inset_values.min(), inset_values.max()
axin2 = ax2.inset_axes([0.2, 0.7, 0.6, 0.25], xlim=(x1,x2), ylim=(y1,y2), xticks=[], yticks=[])
ax2.indicate_inset_zoom(axin2, edgecolor="black")
axin2.plot(dpos2, weight)
ax2.set_xlabel("Position (kb)")
ax.set_ylabel("Coverage")
ax2.set_ylabel("PWM score")
fig.savefig("results/hexamer_correlation/hexamer_correlation.png")

fig, (ax1,ax2) = pyplot.subplots(nrows=2)
ax1.plot(pos, small_variations, color="grey")
ax2.plot(dpos2, np.cumsum(weight - weight.mean()), color="green")
fig.savefig("results/hexamer_correlation/hexamer_correlation.png")

x1,x2 = 0.2,0.4
fig, ax = pyplot.subplots(constrained_layout=True)
ax.plot(dpos1, diff, color="black")
ax.set_xlim(x1,x2)
ax.set_xlabel("Position (kb)")
ax2 = ax.twinx()
ax2.plot(dpos2, weight, color="red", alpha=0.7)
ax.set_ylabel("Coverage difference")
ax2.set_ylabel("PWM Score")
ax2.tick_params(axis='y', colors="red")
fig.savefig("results/hexamer_correlation/hexamer_correlation.png")
