UHR = pl.read_csv(
    "results/transcript_coverage/UHR.transcript_coverage.txt.gz", separator="\t"
)
HBRR = pl.read_csv(
    "results/transcript_coverage/HBRR.transcript_coverage.txt.gz", separator="\t"
)

transcripts_in_both = (
    UHR.select("gene")
    .unique()
    .join(HBRR.select("gene").unique(), how="inner", on="gene")
)

UHR = UHR.join(transcripts_in_both, how="inner", on="gene")
HBRR = HBRR.join(transcripts_in_both, how="inner", on="gene")

sample_info = pl.concat(
    [
        pl.read_csv("results/UHR.sample_info.txt", separator="\t"),
        pl.read_csv("results/HBRR.sample_info.txt", separator="\t"),
    ]
)


def compute_corr(cov1, cov2):
    # compute correlation of coverage profiles in each gene
    # for the two samples cov1, cov2
    return (
        cov1.join(cov2, on=["gene", "pos"], suffix="_2")
        .group_by("gene")
        .agg(
            corr=pl.corr(
                "cov",
                "cov_2",
            )
        )
    )


# Compute correlations within tissues
# do this for samples within the same site only
temp = []
for tissue, cov in [("UHR", UHR), ("HBRR", HBRR)]:
    for study in sample_info.filter(tissue=tissue)["study"].unique():
        site = study.split("_")[1]
        sample_ids = sample_info.filter(tissue=tissue, study=study)["ID"]
        for sample1 in sample_ids:
            for sample2 in sample_ids:
                if sample1 <= sample2:
                    continue
                temp.append(
                    compute_corr(
                        cov.filter(sample=sample1),
                        cov.filter(sample=sample2),
                    ).with_columns(
                        sample1=pl.lit(sample1),
                        sample2=pl.lit(sample2),
                        site=pl.lit(site),
                        tissue=pl.lit(tissue),
                    )
                )
intra_corr = pl.concat(temp)

# Compute correlations between tissues
# do this for samples within the same site only
temp = []
for site in ["BGI", "MAY", "CNL"]:
    UHR_study_id = f"PRJNA208369_{site}"
    HBRR_study_id = f"PRJNA208369_{site}_B"
    UHR_samples = sample_info.filter(tissue="UHR", study=UHR_study_id)["ID"]
    HBRR_samples = sample_info.filter(tissue="HBRR", study=HBRR_study_id)["ID"]
    for sample1 in UHR_samples:
        for sample2 in HBRR_samples:
            temp.append(
                compute_corr(
                    UHR.filter(sample=sample1),
                    HBRR.filter(sample=sample2),
                ).with_columns(
                    sample1=pl.lit(sample1),
                    sample2=pl.lit(sample2),
                    site=pl.lit(site),
                )
            )
inter_corr = pl.concat(temp)
