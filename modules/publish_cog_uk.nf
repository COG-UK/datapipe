#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
includeConfig 'config/base.config'


process uk_geography {
    /**
    * Cleans up geography
    * @input uk_fasta, uk_metadata
    * @output uk_metadata_updated
    * @params geography_utils
    */

    input:
    file uk_fasta
    file uk_metadata

    output:
    path "geography/geography.csv"

    script:
    """
    mkdir geography
    mkdir geography_tmp

    fastafunk fetch \
      --in-fasta ${uk_fasta} \
      --in-metadata ${uk_metadata} \
      --index-column sequence_name \
      --filter-column central_sample_id sequence_name sample_date epi_week \
                      adm0 adm1 adm2 adm2_private \
      --out-fasta geography_tmp/fetch.fa \
      --out-metadata geography_tmp/fetch.csv \
      --restrict

    $project_dir/../bin/geography_cleaning.py \
      --metadata geography_tmp/fetch.csv \
      --country-col adm0 \
      --adm1-col adm1 \
      --adm2-col adm2 \
      --outer-postcode-col adm2_private \
      --mapping-utils-dir ${geography_utils} \
      --outdir geography

    rm -rf geography_tmp
    """
}


workflow publish_cog_uk {
    take:
        uk_unaligned_fasta
        uk_full_aligned_fasta
        uk_filtered_aligned_fasta
        uk_metadata
    main:

    emit:

}


workflow {
    publish_cog_uk(params.uk_fasta,
                   params.uk_metadata)
}