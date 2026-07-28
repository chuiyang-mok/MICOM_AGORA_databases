nextflow.enable.dsl=2

params.out = "${launchDir}"
params.manifest = "${params.out}/data/agora201.csv"
params.agora_version = "201"
params.refseq = "232"
params.gtdb = "232"
params.database_version = "1"
params.models = "${params.out}/models"

workflow {
    levels = channel.from(["genus", "species"])

    DownloadGtdbTables()

    refseq_manifests = ConvertToRefseq(levels)
    gtdb_manifests = ConvertToGtdb(refseq_manifests, DownloadGtdbTables.out)

    refseq_manifests.concat(gtdb_manifests).view()

    BuildDB(refseq_manifests.concat(gtdb_manifests)) | GetManifest
}


// Process definitions

process ConvertToRefseq {
    cpus 1
    memory 4.GB
    time 1.h

    input:
    val(level)

    output:
    tuple val(level), val("refseq"), val("${params.refseq}"), path("agora2_refseq_${level}.tsv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from os import path

    manifest = pd.read_csv("${params.manifest}")
    del manifest["sbml_model"]
    manifest["file"] = [path.join("${params.models}", id + ".xml.gz") for id in manifest["MicrobeID"]]
    ranks = pd.Index(["Strain", "Species", "Genus" , "Family","Order", "Class", "Phylum", "Kingdom"][::-1])
    manifest.rename(columns=dict(zip(ranks, ranks.str.lower())), inplace=True)
    manifest.rename(columns={"MicrobeID": "id"}, inplace=True)
    print(manifest)
    manifest.to_csv("agora2_refseq_${level}.tsv", index=False, sep="\\t")
    """
}

process DownloadGtdbTables {
    cpus 1
    memory 8.GB
    time 8.h
    errorStrategy 'retry'
    maxRetries 3

    output:
    path("*.tsv")

    script:
    """
    wget --no-check-certificate -O gtdb_bac.tsv.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/bac120_metadata_r${params.gtdb}.tsv.gz
    wget --no-check-certificate -O gtdb_ar.tsv.gz https://data.gtdb.ecogenomic.org/releases/release${params.gtdb}/${params.gtdb}.0/ar53_metadata_r${params.gtdb}.tsv.gz
    gunzip -f gtdb_bac.tsv.gz
    gunzip -f gtdb_ar.tsv.gz
    """
}

process ConvertToGtdb {
    cpus 1
    memory 4.GB
    time 1.h

    input:
    tuple val(level), val(db), val(ver), path(manifest)
    path(tables)

    output:
    tuple val(level), val("gtdb"), val("${params.gtdb}"), path("agora2_gtdb_${level}.tsv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    gtdb_rank_names = pd.Index(["domain", "phylum", "class", "order", "family", "genus", "species"])

    files = "${tables}".split()
    meta = pd.concat(
        pd.read_csv(fi, sep="\\t", usecols=["gtdb_taxonomy", "ncbi_taxid"])
        for fi in files
    )

    rank_idx = gtdb_rank_names.get_loc("${level}")
    tax = meta.gtdb_taxonomy.str.split(";", expand=True)
    tax.columns = gtdb_rank_names
    meta = pd.concat([meta, tax], axis = 1).drop(["gtdb_taxonomy"], axis=1).drop_duplicates()

    agora = pd.read_csv("${manifest}", sep="\\t").drop(
        ["strain", "species", "genus" , "family","order", "class", "phylum", "kingdom"],
        axis=1, errors="ignore"
    )
    merged = pd.merge(agora, meta, left_on="NCBI Taxonomy ID", right_on="ncbi_taxid")
    rank_mappings = merged.groupby("file")["${level}"].nunique()
    valid = rank_mappings.index[rank_mappings == 1]
    merged = merged[merged.file.isin(valid)].drop_duplicates(subset=["file", "${level}"])
    merged.rename(columns={"domain": "kingdom"}, inplace=True)

    print(f"rank: ${level} - matched: {merged.shape[0]}/{agora.shape[0]}")

    merged.to_csv("agora2_gtdb_${level}.tsv", index=False, sep="\\t")
    """

}

process BuildDB {
    publishDir "${params.out}/databases", mode: "copy", overwrite: true
    cpus 12
    memory 16.GB
    time 8.h

    input:
    tuple val(level), val(db), val(ver), path(manifest)

    output:
    tuple val(level), val(db), val(ver), path("*.qza")

    script:
    """
    qiime micom db --m-meta-file ${manifest} \
        --p-rank ${level} \
        --p-threads ${task.cpus} \
        --verbose \
        --o-metabolic-models agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.qza
    """
}

process GetManifest {
    cpus 1
    memory 4.GB
    time 1.h
    publishDir "${params.out}/manifests", mode: "copy", overwrite: true

    input:
    tuple val(level), val(db), val(ver), path(arti)

    output:
    path("agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.tsv")

    script:
    """
    #!/usr/bin/env python
    from micom.qiime_formats import load_qiime_manifest

    manifest = load_qiime_manifest("${arti}")
    manifest.to_csv(
        "agora${params.agora_version}_${db}${ver}_${level}_${params.database_version}.tsv",
        sep="\\t",
        index=False
    )
    """
}
