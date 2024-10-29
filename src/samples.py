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
    # Moreover, they split the lanes up over multiple SRX entries and therefore we have to manually provide all the SRR values to use
    {"ID": "SRX302130", "study": "SRX302209_BGI", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR896663", "SRR896664", "SRR896665", "SRR896666", "SRR896667", "SRR896668", "SRR896669", "SRR896670", "SRR896671", "SRR896672", "SRR896673", "SRR896674", "SRR896675", "SRR896676", "SRR896677", "SRR896678"]},
    {"ID": "SRX302146", "study": "SRX302209_BGI", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR896679", "SRR896680", "SRR896681", "SRR896682", "SRR896683", "SRR896684", "SRR896685", "SRR896686", "SRR896687", "SRR896688", "SRR896689", "SRR896690", "SRR896691", "SRR896692", "SRR896693", "SRR896694"]},
    {"ID": "SRX302162", "study": "SRX302209_BGI", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR896695", "SRR896696", "SRR896697", "SRR896698", "SRR896699", "SRR896700", "SRR896701", "SRR896702", "SRR896703", "SRR896704", "SRR896705", "SRR896706", "SRR896707", "SRR896708", "SRR896709", "SRR896710"]},
    {"ID": "SRX302178", "study": "SRX302209_BGI", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR896711", "SRR896712", "SRR896713", "SRR896714", "SRR896715", "SRR896716", "SRR896717", "SRR896718", "SRR896719", "SRR896720", "SRR896721", "SRR896722", "SRR896723", "SRR896724", "SRR896725", "SRR896726"]},
    {"ID": "SRX302194", "study": "SRX302209_BGI", "library_id": 5, "tissue": "UHR", "all_SRR": ["SRR896727", "SRR896728", "SRR896729", "SRR896730", "SRR896731", "SRR896732", "SRR896733", "SRR896734", "SRR896735", "SRR896736", "SRR896737", "SRR896738", "SRR896739", "SRR896740", "SRR896741", "SRR896742"]},

    {"ID": "SRX302514", "study": "SRX302209_CNL", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR897047", "SRR897048", "SRR897049", "SRR897050", "SRR897051", "SRR897052", "SRR897053", "SRR897054", "SRR897055", "SRR897056", "SRR897057", "SRR897058", "SRR897059", "SRR897060", "SRR897061"]},
    {"ID": "SRX302529", "study": "SRX302209_CNL", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR897062", "SRR897063", "SRR897064", "SRR897065", "SRR897066", "SRR897067", "SRR897068", "SRR897069", "SRR897070", "SRR897071", "SRR897072", "SRR897073", "SRR897074", "SRR897075", "SRR897076"]},
    {"ID": "SRX302544", "study": "SRX302209_CNL", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR897077", "SRR897078", "SRR897079", "SRR897080", "SRR897081", "SRR897082", "SRR897083", "SRR897084", "SRR897085", "SRR897086", "SRR897087", "SRR897088", "SRR897089", "SRR897090", "SRR897091"]},
    {"ID": "SRX302559", "study": "SRX302209_CNL", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR897092", "SRR897093", "SRR897094", "SRR897095", "SRR897096", "SRR897097", "SRR897098", "SRR897099", "SRR897100", "SRR897101", "SRR897102", "SRR897103", "SRR897104", "SRR897105", "SRR897106"]},
    {"ID": "SRX302574", "study": "SRX302209_CNL", "library_id": 5, "tissue": "UHR", "all_SRR": ["SRR897107", "SRR897108", "SRR897109", "SRR897110", "SRR897111", "SRR897112", "SRR897113", "SRR897114", "SRR897115", "SRR897116", "SRR897117", "SRR897118", "SRR897119", "SRR897120", "SRR897121"]},

    {"ID": "SRX302874", "study": "SRX302209_MAY", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR897407", "SRR897408", "SRR897409", "SRR897410", "SRR897411", "SRR897412", "SRR897413", "SRR897414", "SRR897415", "SRR897416", "SRR897417", "SRR897418", "SRR897419", "SRR897420", "SRR897421", "SRR897422"]},
    {"ID": "SRX302890", "study": "SRX302209_MAY", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR897423", "SRR897424", "SRR897425", "SRR897426", "SRR897427", "SRR897428", "SRR897429", "SRR897430", "SRR897431", "SRR897432", "SRR897433", "SRR897434", "SRR897435", "SRR897436", "SRR897437", "SRR897438"]},
    {"ID": "SRX302906", "study": "SRX302209_MAY", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR897439", "SRR897440", "SRR897441", "SRR897442", "SRR897443", "SRR897444", "SRR897445", "SRR897446", "SRR897447", "SRR897448", "SRR897449", "SRR897450", "SRR897451", "SRR897452", "SRR897453", "SRR897454"]},
    {"ID": "SRX302922", "study": "SRX302209_MAY", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR897455", "SRR897456", "SRR897457", "SRR897458", "SRR897459", "SRR897460", "SRR897461", "SRR897462", "SRR897463", "SRR897464", "SRR897465", "SRR897466", "SRR897467", "SRR897468", "SRR897469", "SRR897470"]},
    {"ID": "SRX302938", "study": "SRX302209_MAY", "library_id": 5, "tissue": "UHR", "all_SRR": ["SRR897471", "SRR897472", "SRR897473", "SRR897474", "SRR897475", "SRR897476", "SRR897477", "SRR897478", "SRR897479", "SRR897480", "SRR897481", "SRR897482", "SRR897483", "SRR897484", "SRR897485", "SRR897486"]},

    #{"ID": "SRX303258", "study": "SRX302209_NVS", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR897791", "SRR897792", "SRR897793", "SRR897794", "SRR897795", "SRR897796", "SRR897797", "SRR897798", "SRR897799", "SRR897800", "SRR897801", "SRR897802", "SRR897803", "SRR897804", "SRR897805", "SRR897806"]},
    #{"ID": "SRX303274", "study": "SRX302209_NVS", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR897807", "SRR897808", "SRR897809", "SRR897810", "SRR897811", "SRR897812", "SRR897813", "SRR897814", "SRR897815", "SRR897816", "SRR897817", "SRR897818", "SRR897819", "SRR897820", "SRR897821", "SRR897822"]},
    #{"ID": "SRX303290", "study": "SRX302209_NVS", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR897823", "SRR897824", "SRR897825", "SRR897826", "SRR897827", "SRR897828", "SRR897829", "SRR897830", "SRR897831", "SRR897832", "SRR897833", "SRR897834", "SRR897835", "SRR897836", "SRR897837", "SRR897838"]},
    #{"ID": "SRX303306", "study": "SRX302209_NVS", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR897839", "SRR897840", "SRR897841", "SRR897842", "SRR897843", "SRR897844", "SRR897845", "SRR897846", "SRR897847", "SRR897848", "SRR897849", "SRR897850", "SRR897851", "SRR897852", "SRR897853", "SRR897854"]},
    #{"ID": "SRX510226", "study": "SRX302209_NYG", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR1215996", "SRR1215997", "SRR1215998", "SRR1215999", "SRR1216000", "SRR1216001", "SRR1216002"]},
    #{"ID": "SRX510233", "study": "SRX302209_NYG", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR1216003", "SRR1216004", "SRR1216005", "SRR1216006", "SRR1216007", "SRR1216008", "SRR1216009"]},
    #{"ID": "SRX510240", "study": "SRX302209_NYG", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR1216010", "SRR1216011", "SRR1216012", "SRR1216013", "SRR1216014", "SRR1216015", "SRR1216016"]},
    #{"ID": "SRX510247", "study": "SRX302209_NYG", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR1216017", "SRR1216018", "SRR1216019", "SRR1216020", "SRR1216021", "SRR1216022", "SRR1216023"]}
    #{"ID": "SRX303706", "study": "SRX302209_AGR", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR898239", "SRR898240", "SRR898241", "SRR898242", "SRR898243", "SRR898244", "SRR898245", "SRR898246", "SRR898247", "SRR898248", "SRR898249", "SRR898250", "SRR898251", "SRR898252", "SRR898253", "SRR898254"]},
    #{"ID": "SRX303722", "study": "SRX302209_AGR", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR898255", "SRR898256", "SRR898257", "SRR898258", "SRR898259", "SRR898260", "SRR898261", "SRR898262", "SRR898263", "SRR898264", "SRR898265", "SRR898266", "SRR898267", "SRR898268", "SRR898269", "SRR898270"]},
    #{"ID": "SRX303738", "study": "SRX302209_AGR", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR898271", "SRR898272", "SRR898273", "SRR898274", "SRR898275", "SRR898276", "SRR898277", "SRR898278", "SRR898279", "SRR898280", "SRR898281", "SRR898282", "SRR898283", "SRR898284", "SRR898285", "SRR898286"]},
    #{"ID": "SRX303754", "study": "SRX302209_AGR", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR898287", "SRR898288", "SRR898289", "SRR898290", "SRR898291", "SRR898292", "SRR898293", "SRR898294", "SRR898295", "SRR898296", "SRR898297", "SRR898298", "SRR898299", "SRR898300", "SRR898301", "SRR898302"]},
    #{"ID": "SRX303578", "study": "SRX302209_COH", "library_id": 1, "tissue": "UHR", "all_SRR": ["SRR898111", "SRR898112", "SRR898113", "SRR898114", "SRR898115", "SRR898116", "SRR898117", "SRR898118"]},
    #{"ID": "SRX303586", "study": "SRX302209_COH", "library_id": 2, "tissue": "UHR", "all_SRR": ["SRR898119", "SRR898120", "SRR898121", "SRR898122", "SRR898123", "SRR898124", "SRR898125", "SRR898126"]},
    #{"ID": "SRX303594", "study": "SRX302209_COH", "library_id": 3, "tissue": "UHR", "all_SRR": ["SRR898127", "SRR898128", "SRR898129", "SRR898130", "SRR898131", "SRR898132", "SRR898133", "SRR898134"]},
    #{"ID": "SRX303602", "study": "SRX302209_COH", "library_id": 4, "tissue": "UHR", "all_SRR": ["SRR898135", "SRR898136", "SRR898137", "SRR898138", "SRR898139", "SRR898140", "SRR898141", "SRR898142"]},
]

for sample in samples:
    if "PCR_cycle_count" not in sample:
        sample["PCR_cycle_count"] = "unknown"
    sample['name'] = f"{sample['ID']}.PCR={sample['PCR_cycle_count']}"

# Check all unique
assert len(set([s['ID'] for s in samples])) == len(samples)
