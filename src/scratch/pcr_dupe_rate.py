import pathlib
import sqlite3
import scipy.ndimage
import scipy.stats
import numpy as np
import polars as pl

library_info = (
    pl.read_excel(
        "../SAMPLE_METADATA/Grant_lab_RNA-seq_libraries_summary_20211101.xlsx"
    )
    .with_columns(library_num=pl.col("LABEL").str.strip_prefix("Lib"))
    .with_columns(
        library_id=pl.col("LABEL")
        + "-"
        + pl.col("SUBLABEL").str.replace("_", "-")
        + "_S"
        + pl.col("library_num"),
        selection_method=pl.col("Selection").replace(
            {
                "PolyA selection": "PolyA",
                "Old Ribo Depletion": "OldRibo",
                "Nu Ribo Depletion": "NuRibo",
                "No selection": "NoSelect",
            }
        ),
    )
)
library_ids = {
    "NoSelect": library_info.filter(selection_method="NoSelect")["library_id"],
    "PolyA": library_info.filter(selection_method="PolyA")["library_id"],
    "Ribo pull-down": library_info.filter(selection_method="OldRibo")["library_id"],
    "Ribo digestion": library_info.filter(selection_method="NuRibo")["library_id"],
}

all_data = pl.concat(

temp = []
for (selection_method, transcript_id), transcript_data in summed_data.group_by(
    ["selection_method", "transcript_id"]
):
    cov = transcript_data["cov"].to_numpy()
    dupe_rate = transcript_data["dupe_rate"].to_numpy()
    mask = np.array(transcript_data["mask"]).astype(bool)

    either_mask = mask[:-1] | mask[1:]

    corr = np.corrcoef(cov[mask], dupe_rate[mask])[0, 1]
    diff_corr = np.corrcoef(np.diff(cov)[either_mask], np.diff(dupe_rate)[either_mask])[
        0, 1
    ]

    smoothed_cov = scipy.ndimage.gaussian_filter1d(cov, sigma=30, mode="constant")
    smoothed_dupe_rate = scipy.ndimage.gaussian_filter1d(
        dupe_rate,
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
            "selection_method": selection_method,
            "transcript_id": transcript_id,
            "corr": corr,
            "diff_corr": diff_corr,
            "smoothed_corr": smoothed_corr,
            "smoothed_diff_corr": smoothed_diff_corr,
        }
    )

results = pl.DataFrame(temp).melt(
    id_vars=["selection_method", "transcript_id"],
    variable_name="corr type",
    value_name="corr",
)


def wilcox(data):
    return scipy.stats.wilcoxon(data).pvalue


print(
    results.group_by(["selection_method", "corr type"]).agg(
        pl.col("corr").median().alias("median corr"),
        pl.col("corr").min().alias("min corr"),
        pl.col("corr").max().alias("max corr"),
        pl.col("corr").map_batches(wilcox).alias("wilcox p"),
    )
)

from matplotlib import pyplot

fig, ax = pyplot.subplots(figsize=(4, 2), sharey=True, constrained_layout=True)
for i, selection_method in enumerate(library_ids.keys()):
    data = results.filter(
        pl.col("corr type") == "smoothed_diff_corr",
        selection_method=selection_method,
    )
    ax.scatter(i + np.random.random(size=len(data["corr"])) * 0.6 - 0.3, data["corr"])
ax.set_xticks(range(len(library_ids)))
ax.set_xticklabels(library_ids.keys())
ax.set_ylim(-1, 1)
ax.set_ylabel("correlation")
fig.savefig("../RESULTS/scratch/pcr_dupe_rate.png", dpi=300)
