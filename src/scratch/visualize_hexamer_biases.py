import polars as pl

samples = [
 {'sample_id': 'SRX16386863', 'study': 'GSE208768'},
 {'sample_id': 'SRX4080514', 'study': 'GSE77221'},
 {'sample_id': 'SRX16386864', 'study': 'GSE208768'},
 {'sample_id': 'SRX4080520', 'study': 'GSE77221'},
 {'sample_id': 'SRX13396189', 'study': 'PRJNA788430'},
 {'sample_id': 'SRX14468347', 'study': 'PRJNA816471'},
 {'sample_id': 'SRX4393369', 'study': 'GSE117134'},
 {'sample_id': 'SRX13396186', 'study': 'PRJNA788430'},
 {'sample_id': 'SRX14468350', 'study': 'PRJNA816471'},
 {'sample_id': 'SRX4393368', 'study': 'GSE117134'},
]
studies = {s['sample_id']: s['study'] for s in samples}

sample_ids = [s['sample_id'] for s in samples]

all_hexamers = []
for sample_id in sample_ids:
    study = studies[sample_id]
    hex = pl.read_csv(f"data/{sample_id}/hexamer_bias.txt", separator="\t").with_columns(sample_id = pl.lit(sample_id), study = pl.lit(study))
    all_hexamers.append(hex)
hexamers = pl.concat(all_hexamers) \
    .with_columns(
        tot = pl.col('A') + pl.col('C') + pl.col('G') + pl.col('T')
    ).with_columns(
        A_freq = pl.col('A') / pl.col('tot'),
        C_freq = pl.col('C') / pl.col('tot'),
        G_freq = pl.col('G') / pl.col('tot'),
        T_freq = pl.col('T') / pl.col('tot'),
    )
