import pathlib
import sqlite3
import numpy as np
import polars as pl
import scipy.ndimage
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as pyplot

cycle_counts_by_id = {
    "SRX3357955": 8,
    "SRX3357956": 9,
    "SRX3357957": 10,
    "SRX3357958": 11,
    "SRX3357959": 13,
}

select_genes = [
    "ENSMUSG00000000594",
    "ENSMUSG00000007836",
    "ENSMUSG00000032437",
]

outdir = pathlib.Path("results/testis/pcr_dupe_rate_coverage_plot/")
outdir.mkdir(exist_ok=True)

data = pl.read_csv(
    "results/testis/pcr_dupe_rate_coverage.aggregated.txt", separator="\t"
)
sample_ids = sorted(data["sample_id"].unique())

annotation_db = "data/Mus_musculus.GRCm38.102.gtf.sqlite"
with sqlite3.connect(annotation_db) as conn:
    # exons = pl.read_database("SELECT * FROM exon", conn)
    # tx2exon = pl.read_database("SELECT * FROM tx2exon", conn)
    gene = pl.read_database("SELECT * FROM gene", conn)
    tx = pl.read_database("SELECT * FROM tx", conn)


###### PLOT SELECT GENES TOGETHER
color_by_cycle_count = {
    cycles: mpl.colormaps["viridis"]((cycles - 8) / 5) for cycles in [8, 9, 10, 11, 13]
}
handles = {}
for format in ["norm", "raw"]:
    fig, axes = pyplot.subplots(
        figsize=(10, 4),
        nrows=2,
        ncols=len(select_genes),
        sharex="col",
        gridspec_kw=dict(right=0.85),
    )
    for i, gene_id in enumerate(select_genes):
        gene_name = gene.filter(gene_id=gene_id)["gene_name"][0]
        cov_ax = axes[0, i]
        dupe_ax = axes[1, i]
        all_cov_and_dupes = data.filter(gene_id=gene_id).with_columns(
            norm_cov=pl.col("cov")
            / pl.col("cov").mean().over(["sample_id", "cycle_count"]),
            norm_dupe_rate=pl.col("dupe_rate")
            / pl.col("dupe_rate").mean().over(["sample_id", "cycle_count"]),
        )

        for (cycle_count,), cov_and_dupe in all_cov_and_dupes.group_by("cycle_count"):
            if format == "raw":
                (h,) = cov_ax.plot(
                    cov_and_dupe["loc"],
                    cov_and_dupe["cov"],
                    color=color_by_cycle_count[cycle_count],
                    label="coverage",
                )

                mask = np.array(cov_and_dupe["mask"]).astype(bool)
                masked_dupe_rate = np.array(cov_and_dupe["dupe_rate"]).astype(float)
                masked_dupe_rate[~mask] = float("nan")
                dupe_ax.plot(
                    cov_and_dupe["loc"],
                    masked_dupe_rate,
                    color=color_by_cycle_count[cycle_count],
                    label="dupe rate",
                )
            else:
                (h,) = cov_ax.plot(
                    cov_and_dupe["loc"],
                    cov_and_dupe["norm_cov"],
                    color=color_by_cycle_count[cycle_count],
                    label="coverage",
                )

                mask = np.array(cov_and_dupe["mask"]).astype(bool)
                masked_dupe_rate = np.array(cov_and_dupe["norm_dupe_rate"]).astype(
                    float
                )
                masked_dupe_rate[~mask] = float("nan")
                dupe_ax.plot(
                    cov_and_dupe["loc"],
                    masked_dupe_rate,
                    color=color_by_cycle_count[cycle_count],
                    label="dupe rate",
                )
            handles[cycle_count] = h

        if i == 0:
            cov_ax.set_ylabel("Coverage")
            dupe_ax.set_ylabel("PCR dupe rate")
        cov_ax.set_title(f"{gene_id}\n{gene_name}")
        dupe_ax.set_xlabel("Position")
    fig.legend(
        handles=[handles[c] for c in color_by_cycle_count.keys()],
        labels=[str(c) for c in color_by_cycle_count.keys()],
        loc="center right",
        title="PCR cycle count",
    )
    fig.tight_layout()
    fig.savefig(outdir / f"select_genes.{format}.svg")
    fig.savefig(outdir / f"select_genes.{format}.png")
    pyplot.close()

####### PLOT INDIVIDUAL GENES
# for gene_id, tx_id in transcripts.iter_rows():
#    gene_name = gene.filter(gene_id=gene_id)["gene_name"][0]
#
#    # Make the plot
#    fig, axes = pyplot.subplots(
#        figsize=(5, 1.5 * len(sample_ids)),
#        nrows=len(sample_ids),
#        layout="constrained",
#        sharex=True,
#        sharey=False,
#    )
#
#    for ax, sample_id in zip(axes, sample_ids):
#        both = get_cov_and_dupe_rate(gene_id, sample_id)
#        mask = np.array(both["mask"]).astype(bool)
#        masked_dupe_rate = np.array(both["dupe_rate"]).astype(float)
#        masked_dupe_rate[~mask] = float("nan")
#        (dupe_rate_lines,) = ax.plot(
#            both["loc"], masked_dupe_rate, color="b", label="dupe rate"
#        )
#        ax.set_ylabel(f"{sample_id}\nPCR dupe rate")
#        ax.set_xlabel("loc")
#        ax_true = ax.twinx()
#        (coverage_lines,) = ax_true.plot(
#            both["loc"], both["cov"], color="k", label="coverage"
#        )
#        ax_true.set_ylabel("coverage")
#    fig.legend(handles=[dupe_rate_lines, coverage_lines], loc="upper right")
#    fig.suptitle(f"{gene_id} - {gene_name}")
#    fig.savefig(outdir / f"{gene_id}.{gene_name}.png", dpi=400)
#    pyplot.close()


### COMPUTE CORRELATION STATS
def weighted_gaussian_filter(x, weights, sigma, mode):
    # also pulls towards the mean very slightly
    # we only use the values where cov > 100 so this makes little difference
    pseudo_count = 0.1
    prior = np.sum(x * weights) / np.sum(weights)
    weighted_smoothed = scipy.ndimage.gaussian_filter1d(
        x * weights, sigma=sigma, mode=mode
    )
    weights = scipy.ndimage.gaussian_filter1d(weights, sigma=sigma, mode=mode)
    return (weighted_smoothed + prior * pseudo_count) / (weights + pseudo_count)


temp = []
for (sample_id, transcript_id), transcript_data in data.group_by(
    ["sample_id", "transcript_id"]
):
    cov = transcript_data["cov"].to_numpy().astype(float)
    dupe_rate = transcript_data["dupe_rate"].to_numpy().astype(float)
    mask = np.array(transcript_data["mask"]).astype(bool)

    either_mask = mask[:-1] | mask[1:]

    corr = np.corrcoef(cov[mask], dupe_rate[mask])[0, 1]
    diff_corr = np.corrcoef(
        np.diff(cov)[either_mask],
        np.diff(dupe_rate)[either_mask],
    )[0, 1]

    smoothed_cov = scipy.ndimage.gaussian_filter1d(cov, sigma=30, mode="constant")
    dupe_rate2 = dupe_rate.copy()
    dupe_rate2[np.isnan(dupe_rate2)] = 0
    smoothed_dupe_rate = weighted_gaussian_filter(
        dupe_rate2,
        weights=cov,
        sigma=30,
        mode="constant",
    )
    smoothed_corr = np.corrcoef(smoothed_cov[mask], smoothed_dupe_rate[mask])[0, 1]
    smoothed_diff_corr = np.corrcoef(
        np.diff(smoothed_cov)[either_mask],
        np.diff(smoothed_dupe_rate)[either_mask],
    )[0, 1]
    temp.append(
        {
            "sample_id": sample_id,
            "cycle_count": cycle_counts_by_id[sample_id],
            "transcript_id": transcript_id,
            "corr": corr,
            "diff_corr": diff_corr,
            "smoothed_corr": smoothed_corr,
            "smoothed_diff_corr": smoothed_diff_corr,
        }
    )

results = pl.DataFrame(temp).melt(
    id_vars=["sample_id", "cycle_count", "transcript_id"],
    variable_name="corr type",
    value_name="corr",
)


# def wilcox(data):
#    return scipy.stats.wilcoxon(data).pvalue
#
#
# print(
#    results.group_by(["cycle_count", "corr type"]).agg(
#        pl.col("corr").median().alias("median corr"),
#        pl.col("corr").min().alias("min corr"),
#        pl.col("corr").max().alias("max corr"),
#        pl.col("corr").map_batches(wilcox).alias("wilcox p"),
#    )
# )

### PLOT CORRELATION STATS
fig, ax = pyplot.subplots(figsize=(3, 2), sharey=True, constrained_layout=True)
for i, sample_id in enumerate(sample_ids):
    df = results.filter(
        pl.col("corr type") == "smoothed_diff_corr",
        sample_id=sample_id,
    )
    ax.scatter(i + np.random.random(size=len(df["corr"])) * 0.6 - 0.3, df["corr"])
ax.set_xticks(range(len(sample_ids)))
ax.set_xticklabels([cycle_counts_by_id[sample_id] for sample_id in sample_ids])
ax.set_ylim(-1, 1)
ax.set_ylabel("correlation")
ax.set_xlabel("PCR cycle number")
fig.savefig(outdir / "pcr_dupe_rate_correlation.png", dpi=300)
