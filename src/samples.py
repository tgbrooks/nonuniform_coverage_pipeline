samples = [
    # Bulk Mouse Testis samples with varying PCR cycle counts and input quantity
    {
        "ID": "SRX3357955",
        "PCR_cycle_count": 8,
        "study": "PRJNA416930",
        "tissue": "testis",
    },
    {
        "ID": "SRX3357956",
        "PCR_cycle_count": 9,
        "study": "PRJNA416930",
        "tissue": "testis",
    },
    {
        "ID": "SRX3357957",
        "PCR_cycle_count": 10,
        "study": "PRJNA416930",
        "tissue": "testis",
    },
    {
        "ID": "SRX3357958",
        "PCR_cycle_count": 11,
        "study": "PRJNA416930",
        "tissue": "testis",
    },
    {
        "ID": "SRX3357959",
        "PCR_cycle_count": 13,
        "study": "PRJNA416930",
        "tissue": "testis",
    },

    # Ribo-Zero Mouse Liver studies
    {
        "ID": "SRX4080514",
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
        "ID": "SRX3304756",
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
        "ID": "SRX16386863",
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
        "ID": "SRX13396189",
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
]

for sample in samples:
    sample['name'] = f"{sample['ID']}.PCR={sample['PCR_cycle_count']}"
