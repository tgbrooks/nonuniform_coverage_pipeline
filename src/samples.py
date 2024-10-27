samples = [
    # Bulk Mouse Testis samples with varying PCR cycle counts and input quantity
    {
        "ID": "SRX3357955",
        "PCR_cycle_count": 8,
        "study": "PRJNA416930",
        "tissue": "testis",
        "UMI": { "read1": "NNNNNNNN", "read2": "NNNNNNNN" },
    },
    {
        "ID": "SRX3357956",
        "PCR_cycle_count": 9,
        "study": "PRJNA416930",
        "tissue": "testis",
        "UMI": { "read1": "NNNNNNNN", "read2": "NNNNNNNN" },
    },
    {
        "ID": "SRX3357957",
        "PCR_cycle_count": 10,
        "study": "PRJNA416930",
        "tissue": "testis",
        "UMI": { "read1": "NNNNNNNN", "read2": "NNNNNNNN" },
    },
    {
        "ID": "SRX3357958",
        "PCR_cycle_count": 11,
        "study": "PRJNA416930",
        "tissue": "testis",
        "UMI": { "read1": "NNNNNNNN", "read2": "NNNNNNNN" },
    },
    {
        "ID": "SRX3357959",
        "PCR_cycle_count": 13,
        "study": "PRJNA416930",
        "tissue": "testis",
        "UMI": { "read1": "NNNNNNNN", "read2": "NNNNNNNN" },
    },

    # Ribo-Zero Mouse Liver studies
    {
        "ID": "SRX4080514",
        "study": "GSE77221",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX4080520",
        "study": "GSE77221",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX6902444",
        "study": "GSE138019",
        "PCR_cycle_count": "unknown",
        # Single-end, can't alpine
        "tissue": "liver",
    },
    {
        "ID": "SRX6902450",
        "study": "GSE138019",
        "PCR_cycle_count": "unknown",
        # Single-end, can't alpine
        "tissue": "liver",
    },
    {
        "ID": "SRX3304756",
        "study": "GSE105413",
        "PCR_cycle_count": "unknown",
        # Single-end, can't alpine
        "tissue": "liver",
    },
    {
        "ID": "SRX3304764",
        "study": "GSE105413",
        "PCR_cycle_count": "unknown",
        # Single-end, can't alpine
        "tissue": "liver",
    },
    {
        "ID": "SRX4393368",
        "study": "GSE117134",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX4393369",
        "study": "GSE117134",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX16386863",
        "study": "GSE208768",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX16386864",
        "study": "GSE208768",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX14468350",
        "study": "PRJNA816471",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX14468347",
        "study": "PRJNA816471",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX13396189",
        "study": "PRJNA788430",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX13396186",
        "study": "PRJNA788430",
        "PCR_cycle_count": "unknown",
        "tissue": "liver",
    },
    {
        "ID": "SRX11694510",
        "study": "PRJNA753198",
        "PCR_cycle_count": "unknown",
        "tissue": "liver"
    },
    {
        "ID": "SRX11694499",
        "study": "PRJNA753198",
        "PCR_cycle_count": "unknown",
        "tissue": "liver"
    },

    ## Small RNA samples:
    #{
    #    "ID": "SRX3357954",
    #    "PCR_cycle_count": 30,
    #},
    #{
    #    "ID": "SRX3357953",
    #    "PCR_cycle_count": 14,
    #},
    #{
    #    "ID": "SRX3357952",
    #    "PCR_cycle_count": 16,
    #},
    #{
    #    "ID": "SRX3357951",
    #    "PCR_cycle_count": 18,
    #},
    #{
    #    "ID": "SRX3357950",
    #    "PCR_cycle_count": 20,
    #},
    #{
    #    "ID": "SRX3357949",
    #    "PCR_cycle_count": 22,
    #},
    #{
    #    "ID": "SRX3357948",
    #    "PCR_cycle_count": 24,
    #},
    #{
    #    "ID": "SRX3357947",
    #    "PCR_cycle_count": 26,
    #},
    #{
    #    "ID": "SRX3357946",
    #    "PCR_cycle_count": 28,
    #},

    ## SEQC samples
    # Library ID of 1,2,3,4 indicates that the libraries were prepared by the central SEQC group
    # and these were sequenced at 6 sites.Library ID of 5 indicates it was a site-specific library.
    # Only the 3 "official" sites prepared the 5th libraries
    {"ID": "SRX302132", "study": "PRJNA208369_BGI", "library_id": 1, "tissue": "UHR"},
    {"ID": "SRX302148", "study": "PRJNA208369_BGI", "library_id": 2, "tissue": "UHR"},
    {"ID": "SRX302164", "study": "PRJNA208369_BGI", "library_id": 3, "tissue": "UHR"},
    {"ID": "SRX302180", "study": "PRJNA208369_BGI", "library_id": 4, "tissue": "UHR"},
    {"ID": "SRX302196", "study": "PRJNA208369_BGI", "library_id": 5, "tissue": "UHR"},

    {"ID": "SRX302516", "study": "PRJNA208369_CNL", "library_id": 1, "tissue": "UHR"},
    {"ID": "SRX302531", "study": "PRJNA208369_CNL", "library_id": 2, "tissue": "UHR"},
    {"ID": "SRX302546", "study": "PRJNA208369_CNL", "library_id": 3, "tissue": "UHR"},
    {"ID": "SRX302561", "study": "PRJNA208369_CNL", "library_id": 4, "tissue": "UHR"},
    {"ID": "SRX302576", "study": "PRJNA208369_CNL", "library_id": 5, "tissue": "UHR"},

    {"ID": "SRX302876", "study": "PRJNA208369_MAY", "library_id": 1, "tissue": "UHR"},
    {"ID": "SRX302892", "study": "PRJNA208369_MAY", "library_id": 2, "tissue": "UHR"},
    {"ID": "SRX302908", "study": "PRJNA208369_MAY", "library_id": 3, "tissue": "UHR"},
    {"ID": "SRX302924", "study": "PRJNA208369_MAY", "library_id": 4, "tissue": "UHR"},
    {"ID": "SRX302940", "study": "PRJNA208369_MAY", "library_id": 5, "tissue": "UHR"},

    #{"ID": "SRX510226", "study": "PRJNA208369_NYG", "library_id": 1, "tissue": "UHR"},
    #{"ID": "SRX510233", "study": "PRJNA208369_NYG", "library_id": 2, "tissue": "UHR"},
    #{"ID": "SRX510240", "study": "PRJNA208369_NYG", "library_id": 3, "tissue": "UHR"},
    #{"ID": "SRX510247", "study": "PRJNA208369_NYG", "library_id": 4, "tissue": "UHR"},

    #{"ID": "SRX303260", "study": "PRJNA208369_NVS", "library_id": 1, "tissue": "UHR"},
    #{"ID": "SRX303276", "study": "PRJNA208369_NVS", "library_id": 2, "tissue": "UHR"},
    #{"ID": "SRX303292", "study": "PRJNA208369_NVS", "library_id": 3, "tissue": "UHR"},
    #{"ID": "SRX303308", "study": "PRJNA208369_NVS", "library_id": 4, "tissue": "UHR"},

    #{"ID": "SRX303579", "study": "PRJNA208369_COH", "library_id": 1, "tissue": "UHR"},
    #{"ID": "SRX303587", "study": "PRJNA208369_COH", "library_id": 2, "tissue": "UHR"},
    #{"ID": "SRX303595", "study": "PRJNA208369_COH", "library_id": 3, "tissue": "UHR"},
    #{"ID": "SRX303603", "study": "PRJNA208369_COH", "library_id": 4, "tissue": "UHR"},

    #{"ID": "SRX303708", "study": "PRJNA208369_AGR", "library_id": 1, "tissue": "UHR"},
    #{"ID": "SRX303724", "study": "PRJNA208369_AGR", "library_id": 2, "tissue": "UHR"},
    #{"ID": "SRX303740", "study": "PRJNA208369_AGR", "library_id": 3, "tissue": "UHR"},
    #{"ID": "SRX303756", "study": "PRJNA208369_AGR", "library_id": 4, "tissue": "UHR"},
]

for sample in samples:
    if "PCR_cycle_count" not in sample:
        sample["PCR_cycle_count"] = "unknown"
    sample['name'] = f"{sample['ID']}.PCR={sample['PCR_cycle_count']}"

# Check all unique
assert len(set([s['ID'] for s in samples])) == len(samples)
