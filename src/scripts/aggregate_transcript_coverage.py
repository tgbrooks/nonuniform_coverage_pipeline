import pandas as pd

all = []
for covfile in snakemake.input.cov:
    all.append(pd.read_csv(covfile, sep="\t", index_col=0))
pd.concat(all).to_csv(snakemake.output.cov, sep="\t")
