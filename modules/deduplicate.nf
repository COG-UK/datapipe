#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process annotate_with_unmapped_genome_completeness {
    /**
        * Adds a column to metadata with proportion of genome which is complete
        * @input fasta, metadata
        * @output metadata
        */

        input:
        path fasta
        path metadata

        output:
        path "${metadata.baseName}.annotated.csv"

        script:
        """
        $project_dir/../bin/annotate_with_unmapped_genome_completeness.py \
            --in-fasta ${fasta} \
            --in-metadata ${metadata} \
            --out-metadata "${metadata.baseName}.annotated.csv"

        if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${metadata.baseName}.annotated.csv" | wc -l) ]]
        then
            echo \$(cat "${metadata}" | wc -l)
            echo \$(cat "${metadata.baseName}.annotated.csv" | wc -l)
            exit 1
        fi
        """
}

process uk_remove_duplicates_COGID_by_proportionN {
    /**
    * Where duplicate COGID, keeps the most complete
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.deduplicated_by_cogid.fa", emit: uk_fasta_updated
    path "${uk_metadata.baseName}.deduplicated_by_cogid.csv", emit: uk_metadata_updated

    script:
    """
    $project_dir/../bin/uk_remove_duplicates_COGID_by_proportionN.py \
        --in-fasta ${uk_fasta} \
        --in-metadata ${uk_metadata} \
        --out-fasta "${uk_fasta.baseName}.deduplicated_by_cogid.fa" \
        --out-metadata "${uk_metadata.baseName}.deduplicated_by_cogid.csv"

    if [[ \$(cat "${uk_metadata}" | wc -l) != \$(cat "${uk_metadata.baseName}.deduplicated_by_cogid.csv" | wc -l) ]]
    then
        echo \$(cat "${uk_metadata}" | wc -l)
        echo \$(cat "${uk_metadata.baseName}.deduplicated_by_cogid.csv" | wc -l)
        exit 1
    fi
    """
}


process remove_duplicates_by_date {
    /**
    * Where duplicate sequence_name, keeps the earliest
    * @input fasta, metadata
    * @output fasta_updated, metadata_updated
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.deduplicated.fa", emit: fasta_updated
    path "${metadata.baseName}.deduplicated.csv", emit: metadata_updated

    script:
    """
    $project_dir/../bin/remove_duplicates_by_date.py \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-fasta "${fasta.baseName}.deduplicated.fa" \
        --out-metadata "${metadata.baseName}.deduplicated.csv"

    if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${metadata.baseName}.deduplicated.csv" | wc -l) ]]
    then
        echo \$(cat "${metadata}" | wc -l)
        echo \$(cat "${metadata.baseName}.deduplicated.csv" | wc -l)
        exit 1
    fi
    """
}


process unify_headers {
    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.UH.fa"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${fasta.baseName}.UH.fa", "w") as fasta_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        id_key = "fasta_header"
        if "edin_header" in reader.fieldnames:
            id_key = "edin_header"
        for row in reader:
            if row["why_excluded"]:
                print("excluded")
                continue
            if row[id_key] in alignment:
                record = alignment[row[id_key]]
                fasta_out.write(">" + row["sequence_name"] + "\\n")
                fasta_out.write(str(record.seq) + "\\n")
            else:
                print(id_key, row[id_key])
    """
}


process uk_label_sourceid_duplicates_to_omit {
    /**
    * Where duplicate source_id, labels all but the earliest as duplicates
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    */

    publishDir "${publish_dev}/cog_gisaid/", pattern: "*.log", mode: 'copy'

    input:
    path uk_metadata

    output:
    path "${uk_metadata.baseName}.deduplicated_by_sourceid.csv", emit: uk_metadata_updated
    path "deduplicated_by_sourceid.log", emit: deduplicate_log

    script:
    """
    $project_dir/../bin/uk_label_sourceid_duplicates_to_omit.py \
        --in-metadata ${uk_metadata} \
        --out-metadata "${uk_metadata.baseName}.deduplicated_by_sourceid.csv"

    if [[ \$(cat "${uk_metadata}" | wc -l) != \$(cat "${uk_metadata.baseName}.deduplicated_by_sourceid.csv" | wc -l) ]]
    then
        echo \$(cat "${uk_metadata}" | wc -l)
        echo \$(cat "${uk_metadata.baseName}.deduplicated_by_sourceid.csv" | wc -l)
        exit 1
    fi
    """
}


workflow deduplicate_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        annotate_with_unmapped_genome_completeness(uk_fasta, uk_metadata)
        uk_remove_duplicates_COGID_by_proportionN(uk_fasta, annotate_with_unmapped_genome_completeness.out)
        unify_headers(uk_remove_duplicates_COGID_by_proportionN.out.uk_fasta_updated, uk_remove_duplicates_COGID_by_proportionN.out.uk_metadata_updated)
        uk_label_sourceid_duplicates_to_omit(uk_remove_duplicates_COGID_by_proportionN.out.uk_metadata_updated)
    emit:
        fasta = unify_headers.out
        metadata = uk_label_sourceid_duplicates_to_omit.out.uk_metadata_updated
}


workflow deduplicate_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
    main:
        annotate_with_unmapped_genome_completeness(gisaid_fasta, gisaid_metadata)
        remove_duplicates_by_date(gisaid_fasta, annotate_with_unmapped_genome_completeness.out)
        unify_headers(remove_duplicates_by_date.out.fasta_updated, remove_duplicates_by_date.out.metadata_updated)
    emit:
        fasta = unify_headers.out
        metadata = remove_duplicates_by_date.out.metadata_updated
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    deduplicate_cog_uk(uk_fasta, uk_metadata)
}
