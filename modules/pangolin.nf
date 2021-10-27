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
    if ( params.auto_update_pangolin )
        """
        PANGO_VERSION=\$(pangolin -pv)
        echo \$PANGO_VERSION
        pangolin --update
        sleep 5s
        NEW_PANGO_VERSION=\$(pangolin -pv)
        echo \$NEW_PANGO_VERSION
        if [ "\$PANGO_VERSION" == "\$NEW_PANGO_VERSION" ]; then
            PANGOLIN_UPDATED=false
        else
            PANGOLIN_UPDATED=true
        fi
        """
    else
        """
        PANGOLIN_UPDATED=false
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
            if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${metadata.baseName}.with_previous.csv" | wc -l) ]]
            then
                echo \$(cat "${metadata}" | wc -l)
                echo \$(cat "${metadata.baseName}.with_previous.csv" | wc -l)
                exit 1
            fi
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
    if (params.skip_designation_hash)
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp \
            --skip-designation-hash
        """
    else
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp
        """
}

process run_pangolin_usher {
    /**
    * Runs PANGOLIN on input fasta
    * @input fasta
    * @output pangolin_fasta
    */

    cpus 4

    input:
    path fasta

    output:
    path "pangolin/lineage_report.csv"

    script:
    if (params.skip_designation_hash)
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp \
            --usher \
            -t ${task.cpus} \
            --skip-designation-hash
        """
    else
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp \
            --usher \
            -t ${task.cpus}
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
    path "${metadata.baseName}.with_pangolin.csv", emit: metadata
    path "pango.log", emit: log

    script:
    """
    $project_dir/../bin/prepare_for_pangolin.py \
                  --in-metadata ${metadata} \
                  --previous-metadata ${pangolin_csv} \
                  --out-metadata "${metadata.baseName}.with_pangolin.csv"
    """
}

process add_pangolin_usher_to_metadata {
    /**
    * Adds usher pangolin calls to metadata
    * @input metadata, usher report
    * @output metadata
    */

    input:
    path metadata
    path usher_report

    output:
    path "${metadata.baseName}.with_usher.csv"

    script:
    """
    if [[ \$(head -n1 ${metadata}) == *"fasta_header"* ]]; then
        fastafunk add_columns \
          --in-metadata ${metadata} \
          --in-data ${usher_report} \
          --index-column fasta_header \
          --join-on taxon \
          --new-columns usher_lineage usher_lineages_version \
          --where-column usher_lineage=lineage usher_lineages_version=version \
          --out-metadata "${metadata.baseName}.with_usher.csv"
    else
        fastafunk add_columns \
          --in-metadata ${metadata} \
          --in-data ${usher_report} \
          --index-column edin_header \
          --join-on taxon \
          --new-columns usher_lineage usher_lineages_version \
          --where-column usher_lineage=lineage usher_lineages_version=version \
          --out-metadata "${metadata.baseName}.with_usher.csv"
    fi
    """
}


process announce_summary {
    /**
    * Summarizes pangolin into JSON
    * @input fastas
    */

    input:
    path pango_input
    path pango_log

    output:
    path "announce.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Finished running pangolin ${params.date}*\\n" >> announce.json
                echo "> Number of sequences input to pangolin for new lineage assignments : \$(cat ${pango_input} | grep '>' | wc -l)\\n" >> announce.json
                echo "> \$(cat ${pango_log})\\n" >> announce.json
                echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Finished running pangolin ${params.date}*\\n" >> announce.json
                echo "> Number of sequences input to pangolin for new lineage assignments : \$(cat ${pango_input} | grep '>' | wc -l)\\n" >> announce.json
                echo "> \$(cat ${pango_log})\\n" >> announce.json
                echo '"}' >> announce.json
            """
}

workflow pangolin {
    take:
        in_fasta
        in_metadata
        pangolin_updated
    main:
        extract_sequences_for_pangolin(in_fasta, in_metadata, pangolin_updated)
        extract_sequences_for_pangolin.out.pangolin_fasta.splitFasta( by: params.chunk_size, file: true )
                                                         .set{ pangolin_chunks }
        run_pangolin(pangolin_chunks)
        run_pangolin.out.collectFile(newLine: true, keepHeader: true, skip: 1)
                        .set{ pangolin_result }
        add_new_pangolin_lineages_to_metadata(extract_sequences_for_pangolin.out.metadata_with_previous, pangolin_result)
        if (params.add_usher_pangolin) {
            run_pangolin_usher(pangolin_chunks)
            run_pangolin_usher.out.collectFile(newLine: true, keepHeader: true, skip: 1)
                                  .set{ pangolin_usher_result }
            add_pangolin_usher_to_metadata(add_new_pangolin_lineages_to_metadata.out.metadata, pangolin_usher_result)
            post_pangolin_metadata = add_pangolin_usher_to_metadata.out
        } else {
            post_pangolin_metadata = add_new_pangolin_lineages_to_metadata.out.metadata
        }

        announce_summary(extract_sequences_for_pangolin.out.pangolin_fasta, add_new_pangolin_lineages_to_metadata.out.log)
    emit:
        metadata = post_pangolin_metadata
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    check_for_pangolin_update()

    pangolin(uk_fasta, uk_metadata, check_for_pangolin_update.out)
}
