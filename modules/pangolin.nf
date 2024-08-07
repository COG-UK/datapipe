#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)


process check_for_pangolin_update {
    /**
    * Checks if there is a new version of pangolin and sets param flag if there is
    */
    output:
    env PANGOLIN_UPDATED

    script:
    if ( params.auto_update_pangolin )
        """
        PANGO_VERSION=\$(pangolin --all-versions)
        echo \$PANGO_VERSION
        pangolin --update
        sleep 5s
        NEW_PANGO_VERSION=\$(pangolin --all-versions)
        echo \$NEW_PANGO_VERSION
        if [[ "\$PANGO_VERSION" == "\$NEW_PANGO_VERSION" ]]; then
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

    cpus 4

    input:
    path fasta

    output:
    path "pangolin/lineage_report.csv", emit: report

    script:
    if (params.skip_designation_hash)
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp \
            --analysis-mode pangolearn \
            --skip-designation-hash
        """
    else
        """
        pangolin "${fasta}" \
            --outdir pangolin \
            --tempdir pangolin_tmp \
            -t ${task.cpus} \
            --add-assignment-cache
        """
}


process add_new_pangolin_lineages_to_metadata {
    /**
    * Updates metadata with new PANGOLIN lineage assignments
    * @input metadata, pangolin_csv
    * @output metadata_updated
    */

    memory { task.attempt * metadata.size() * 3.B }

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


process cache_lineages_report {
    /**
    * Creates a map from sequence hash to pangolin report calls
    * @input metadata
    * @output metadata
    */
    publishDir "${publish_dir}/pangolin", pattern: "*.cache.csv", mode: 'copy'

    input:
    path fasta
    path metadata

    output:
    path "${metadata.baseName}.cache.csv", emit: metadata

    script:
    """
    $project_dir/../bin/cache_pangolin_report.py \
        --in-fasta ${fasta} \
        --in-metadata ${metadata} \
        --out-metadata "${metadata.baseName}.cache.csv"
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


process publish_metadata {
    /**
    * Publishes metadata csv for this category
    * @input metadata
    * @output metadata
    */

    publishDir "${publish_dir}", pattern: "*/*.csv", mode: 'copy'

    input:
    path metadata
    val category

    output:
    path "${category}/pangolin_master.csv"

    script:
    """
    mkdir -p ${category}
    cp ${metadata} ${category}/pangolin_master.csv
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
        run_pangolin.out.report.collectFile(newLine: true, keepHeader: true, skip: 1)
                        .set{ pangolin_result }
        post_pangolin_metadata = pangolin_result
        add_new_pangolin_lineages_to_metadata(extract_sequences_for_pangolin.out.metadata_with_previous, post_pangolin_metadata)

        if (params.cache_pangolin){
            cache_lineages_report(in_fasta, post_pangolin_metadata)
        }

        announce_summary(extract_sequences_for_pangolin.out.pangolin_fasta, add_new_pangolin_lineages_to_metadata.out.log)
    emit:
        metadata = add_new_pangolin_lineages_to_metadata.out.metadata
        report = post_pangolin_metadata
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    check_for_pangolin_update()

    pangolin(uk_fasta, uk_metadata, check_for_pangolin_update.out)
    publish_metadata(pangolin.out.report, "pangolin")
}
