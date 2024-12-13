import pathlib
import numpy as np
import polars as pl
from matplotlib import pyplot

data = []

sample_info = pl.read_csv(snakemake.input.sample_info, separator="\t")
study_by_sample = dict(zip(sample_info['ID'], sample_info['study']))

for file in snakemake.input.hexamer:
    sample_id = pathlib.Path(file).parent.name
    study = study_by_sample[sample_id]
    data.append(pl.read_csv(file, separator="\t").with_columns(sample_id=pl.lit(sample_id), study = pl.lit(study)))
hex = pl.concat(data)
hex.write_csv(snakemake.output.tsv, separator="\t")


sample_ids = hex.select('sample_id', 'study').unique().sort('study')['sample_id']
fig, axes = pyplot.subplots(figsize=(10,5), ncols = 5, nrows=2, constrained_layout=True)
for sample_id, ax in zip(sample_ids, axes.T.flatten()):
    study = study_by_sample[sample_id]
    this = hex.filter(sample_id = sample_id).sort("actual_corr")
    x = np.arange(len(this))
    h_ac, = ax.plot(
        x,
        this['actual_corr'],
        color = 'r',
        label="PWM",
    )
    h_rpc = ax.fill_between(
        x,
        this['mean_perm_corr'] - 1.96*this['std_perm_corr'],
        this['mean_perm_corr'] + 1.96*this['std_perm_corr'],
        color = '#888888',
    )
    h_mpc, = ax.plot(
        x,
        this['mean_perm_corr'],
        color = "k",
        label = "Permuted",
    )
    ax.set_title(f"{study}\n{sample_id}")
    if ax in axes[:,0]:
        ax.set_ylabel("correlation")
    else:
        ax.set_yticks([])
    if ax in axes[-1,:]:
        ax.set_xlabel("gene")
fig.legend([h_ac, h_mpc], ["PWM", "permuted"], loc="outside right upper" )
fig.savefig(snakemake.output.pvalue_fig, dpi=400)


print(
    hex.group_by('study', 'sample_id')
    .agg(
        pl.col('actual_corr').median().alias('median corr'),
        pl.col('actual_corr').max().alias('max corr'),
        pl.col('actual_corr').min().alias('min corr'),
        pl.col('p').max().alias('max p')
    )
    .sort('study', 'sample_id')
)
