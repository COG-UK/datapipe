#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process check_for_pangolin_update {
    /**
    * Checks if there is a new version of pangolin and sets param flag if there is
    */
    output:
    env PANGOLIN_UPDATED

    script:
    """
    PANGO_VERSION=\$(pangolin -pv)
    echo \$PANGO_VERSION
    pangolin --update
    NEW_PANGO_VERSION=\$(pangolin -pv)
    echo \$NEW_PANGO_VERSION
    if [ "\$PANGO_VERSION" == "\$NEW_PANGO_VERSION" ]; then
        PANGOLIN_UPDATED=false
    else
        PANGOLIN_UPDATED=true
    fi
    """
}


process extract_sequences_for_pangolin {
    /**
    * If update_all_lineage_assignments flag set, or no previous provided, outputs the input files.
    * Otherwise, extracts lineageless sequences from FASTA to run pangolin on, and updates
    * metadata with previous lineages
    * @input fasta, metadata
    * @output pangolin_fasta, metadata_with_previous
    * @params previous_metadata, update_all_lineage_assignments
    */

    input:
    path fasta
    path metadata
    env PANGOLIN_UPDATED

    output:
    path "${fasta.baseName}.for_pangolin.fa", emit: pangolin_fasta
    path "${metadata.baseName}.with_previous.csv", emit: metadata_with_previous

    script:
    if (params.update_all_lineage_assignments || !params.previous_metadata )
        """
        mv "${fasta}" "${fasta.baseName}.for_pangolin.fa"
        mv "${metadata}" "${metadata.baseName}.with_previous.csv"
        """
    else
        """
        echo "Pangolin updated: \$PANGOLIN_UPDATED"
        if [ \$PANGOLIN_UPDATED == "true" ]
        then
            mv "${fasta}" "${fasta.baseName}.for_pangolin.fa"
            mv "${metadata}" "${metadata.baseName}.with_previous.csv"
        else
            $project_dir/../bin/prepare_for_pangolin.py \
              --in-fasta ${fasta} \
              --in-metadata ${metadata} \
              --previous-metadata ${params.previous_metadata} \
              --out-fasta "${fasta.baseName}.for_pangolin.fa" \
              --out-metadata "${metadata.baseName}.with_previous.csv"
        fi
        """
}

process run_pangolin {
    /**
    * Runs PANGOLIN on input fasta
    * @input fasta
    * @output pangolin_fasta
    */

    input:
    path fasta

    output:
    path "pangolin/lineage_report.csv"

    script:
    """
    pangolin "${fasta}" \
        --outdir pangolin \
        --tempdir pangolin_tmp
    """
}

process add_new_pangolin_lineages_to_metadata {
    /**
    * Updates metadata with new PANGOLIN lineage assignments
    * @input metadata, pangolin_csv
    * @output metadata_updated
    */

    input:
    path metadata
    path pangolin_csv

    output:
    path "${metadata.baseName}.with_pangolin.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    lineage_dict = {}
    with open("${pangolin_csv}", 'r', newline = '') as lineages_in:
        reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")
        for row in reader:
            if row["taxon"] in lineage_dict:
                print("%s occurs more than once in lineages input file" % row["taxon"])
                continue
            lineage_dict[row["taxon"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["pangoLEARN_version"], "probability": row["probability"]}


    with open("${metadata}", 'r', newline = '') as csv_in, \
         open("${metadata.baseName}.with_pangolin.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        column_names = reader.fieldnames + [col for col in ["lineage", "pangoLEARN_version", "probability"] if col not in reader.fieldnames]
        writer = csv.DictWriter(csv_out, fieldnames = column_names, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        if "fasta_header" in reader.fieldnames:
            taxon = "fasta_header"
        else:
            taxon = "edin_header"

        for row in reader:
            fasta_header = row[taxon]
            if fasta_header in lineage_dict:
                row["lineage"] = lineage_dict[fasta_header]["lineage"]
                row["pangoLEARN_version"] = lineage_dict[fasta_header]["pangoLEARN_version"]
                row["probability"] = lineage_dict[fasta_header]["probability"]
            writer.writerow(row)
    """
}

workflow pangolin {
    take:
        in_fasta
        in_metadata
    main:
        check_for_pangolin_update()
        extract_sequences_for_pangolin(in_fasta, in_metadata, check_for_pangolin_update.out)
        run_pangolin(extract_sequences_for_pangolin.out.pangolin_fasta)
        add_new_pangolin_lineages_to_metadata(extract_sequences_for_pangolin.out.metadata_with_previous, run_pangolin.out)
    emit:
        metadata = add_new_pangolin_lineages_to_metadata.out
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)

    pangolin(uk_fasta, uk_metadata)
}
