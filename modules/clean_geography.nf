#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process uk_geography {
    /**
    * Cleans up geography
    * @input uk_fasta, uk_metadata
    * @output geography_metadata
    * @params geography_utils
    */

    memory { 1.GB * task.attempt + uk_fasta.size() * 1.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1

    publishDir "${publish_dev}/", pattern: "geography/*.csv", mode: 'copy'
    publishDir "${publish_dev}/", pattern: "geography/*.txt", mode: 'copy'

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "geography/geography.csv", emit: geography
    path "geography/*.csv"
    path "geography/*.txt"

    script:
    """
    mkdir geography
    mkdir geography_tmp

    fastafunk fetch \
      --in-fasta ${uk_fasta} \
      --in-metadata ${uk_metadata} \
      --index-column sequence_name \
      --filter-column central_sample_id sequence_name sample_date edin_epi_week \
                      adm0 adm1 adm2 adm2_private \
      --out-fasta geography_tmp/fetch.fa \
      --out-metadata geography_tmp/fetch.csv \
      --restrict

    $project_dir/../bin/geography_cleaning/geography_cleaning.py \
      --metadata geography_tmp/fetch.csv \
      --country-col adm0 \
      --adm1-col adm1 \
      --adm2-col adm2 \
      --outer-postcode-col adm2_private \
      --mapping-utils-dir ${geography_utils} \
      --epiweek-col edin_epi_week \
      --outdir geography

    #rm -rf geography_tmp
    """
}


process add_uk_geography_to_metadata {
    /**
    * Adds UK geography to uk metadata
    * @input combined_metadata, geography_metadata
    * @output metadata
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_master.csv"}
    memory { 1.GB * task.attempt + uk_metadata.size() * 2.B }

    input:
    path uk_metadata
    path geography_metadata

    output:
    path "cog_geography.csv", emit: metadata

    script:
    """
    fastafunk add_columns \
          --in-metadata ${uk_metadata} \
          --in-data ${geography_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --force-overwrite \
          --new-columns adm1 adm1_raw adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location safe_location utla utla_code suggested_adm2_grouping \
          --out-metadata "cog_geography.csv"
    """
}


process gisaid_geography {
    /**
    * Cleans up geography
    * @input gisaid_fasta, gisaid_metadata
    * @output geography_metadata
    * @params geography_utils
    */

    memory { 1.GB * task.attempt + fasta.size() * 1.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1

    publishDir "${publish_dev}/", pattern: "geography/*.csv", mode: 'copy'
    publishDir "${publish_dev}/", pattern: "geography/*.txt", mode: 'copy'

    input:
    path gisaid_fasta
    path gisaid_metadata

    output:
    path "geography/geography.csv", emit: geography
    path "geography/*.csv"
    path "geography/*.txt"

    script:
    """
    mkdir geography
    mkdir geography_tmp

    fastafunk fetch \
      --in-fasta ${fasta} \
      --in-metadata ${metadata} \
      --index-column sequence_name \
      --filter-column gisaid_accession sequence_name sample_date epi_week \
                      adm0 adm1 adm2 adm2_private \
      --where-column gisaid_accession=covv_accession_id epi_week=edin_epi_week adm0=edin_admin_0 adm1=edin_admin_1 adm2=edin_admin_2\
      --out-fasta geography_tmp/fetch.fa \
      --out-metadata geography_tmp/fetch.csv \
      --restrict

    $project_dir/../bin/geography_cleaning/geography_cleaning.py \
      --metadata geography_tmp/fetch.csv \
      --country-col adm0 \
      --adm1-col adm1 \
      --adm2-col adm2 \
      --outer-postcode-col adm2_private \
      --mapping-utils-dir ${geography_utils} \
      --epiweek-col epi_week \
      --sample-id-col gisaid_accession \
      --outdir geography

    rm -rf geography_tmp
    """
}


process add_gisaid_geography_to_metadata {
    /**
    * Adds GISAID geography to combined metadata
    * @input gisaid_metadata, geography_metadata
    * @output metadata
    */

    publishDir "${publish_dev}/gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"gisaid_master.csv"}, overwrite: true
    memory { 1.GB * task.attempt + combined_metadata.size() * 2.B }

    input:
    path gisaid_metadata
    path geography_metadata

    output:
    path "gisaid_geography.csv", emit: metadata

    script:
    """
    fastafunk add_columns \
          --in-metadata ${gisaid_metadata} \
          --in-data ${geography_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --force-overwrite \
          --new-columns edin_admin_0 edin_admin_1 edin_admin_2 adm1 adm1_raw adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location safe_location utla utla_code suggested_adm2_grouping \
          --where-column edin_admin_0=adm0 edin_admin_1=adm1 edin_admin_2=adm2 \
          --out-metadata "gisaid_geography.csv"
    """
}


process make_delta_by_utla_summary {
    /**
    * Summarizes delta counts by utla
    * @input metadata
    * @output csv
    */

    publishDir "${publish_dir}/cog", pattern: "*.csv", mode: 'copy', overwrite: false

    input:
    path metadata

    output:
    path "UTLA_genome_counts_${params.date}.csv"

    script:
    """
    $project_dir/../bin/summarise_genomes_by_utla.py \
      --metadata ${metadata} \
      --date ${params.date}
    """
}

process publish_master_metadata {
    /**
    * Publishes master metadata csv for this category
    * @input metadata
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path metadata
    val category

    output:
    path "${category}/${category}_master.csv"

    script:
    """
    mkdir -p ${category}
    cp ${metadata} ${category}/${category}_master.csv
    """
}


geography_utils = file(params.uk_geography)


workflow clean_geography_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        uk_geography(uk_fasta, uk_metadata)
        add_uk_geography_to_metadata(uk_metadata,uk_geography.out.geography)
        make_delta_by_utla_summary(add_uk_geography_to_metadata.out.metadata)
        publish_master_metadata(add_uk_geography_to_metadata.out.metadata, "cog")
    emit:
        metadata = add_uk_geography_to_metadata.out.metadata
}

workflow clean_geography_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
    main:
        if ( params.geography ){
            gisaid_geography(gisaid_fasta, gisaid_metadata)
            add_gisaid_geography_to_metadata(gisaid_metadata,gisaid_geography.out.geography)
            add_gisaid_geography_to_metadata.out.metadata.set{ new_gisaid_metadata }
        } else {
            new_gisaid_metadata = gisaid_metadata
        }
        publish_master_metadata(new_gisaid_metadata, "gisaid")
    emit:
        metadata = new_gisaid_metadata
}

workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    clean_geography_cog_uk(uk_fasta, uk_metadata)
}
