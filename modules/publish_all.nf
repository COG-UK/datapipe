#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process combine_cog_gisaid {
    /**
    * Combines FASTA and METADATA for COG-UK and GISAID
    * @input uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata
    * @output cog_gisaid_fasta, cog_gisaid_metadata
    */

    input:
    path uk_fasta
    path uk_metadata
    path gisaid_fasta
    path gisaid_metadata

    output:
    path "cog_gisaid.fa", emit: fasta
    path "cog_gisaid.csv", emit: metadata

    script:
    """
        fastafunk fetch \
          --in-fasta ${uk_fasta} \
          --in-metadata ${uk_metadata} \
          --index-column sequence_name \
          --filter-column fasta_header covv_accession_id central_sample_id biosample_source_id secondary_identifier root_sample_id \
                          pillar_2 \
                          sequence_name sample_date epi_week \
                          country adm1 adm1_UK adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
                          submission_org_code is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineage_support lineages_version \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target why_excluded \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private lineage_support=probability lineages_version=pangoLEARN_version adm1_UK=adm1\
          --out-fasta "intermediate_cog.fa" \
          --out-metadata "intermediate_cog.csv" \
          --restrict

        fastafunk fetch \
          --in-fasta ${gisaid_fasta} \
          --in-metadata ${gisaid_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id central_sample_id biosample_source_id secondary_identifier root_sample_id \
                          pillar_2 \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
                          submission_org_code is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineage_support lineages_version \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
          --where-column adm1=edin_admin_1 travel_history=edin_travel \
          --out-fasta "intermediate_gisaid.fa" \
          --out-metadata "intermediate_gisaid.csv" \
          --restrict

        fastafunk merge \
          --in-fasta "intermediate_cog.fa" "intermediate_gisaid.fa" \
          --in-metadata "intermediate_cog.csv" "intermediate_gisaid.csv" \
          --out-fasta "cog_gisaid.fa" \
          --out-metadata "cog_gisaid.csv" \
          --index-column sequence_name
    """
}


process combine_variants {
    /**
    * Combines FASTA and METADATA for COG-UK and GISAID
    * @input uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata
    * @output cog_gisaid_fasta, cog_gisaid_metadata
    */

    publishDir "${publish_dev}/COG_GISAID", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_variants.csv"}


    input:
    path uk_fasta
    path uk_variants
    path gisaid_fasta
    path gisaid_variants

    output:
    path "cog_gisaid_variants.csv"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    COG_fasta = SeqIO.index("${uk_fasta}", "fasta")
    GISAID_fasta = SeqIO.index("${gisaid_fasta}", "fasta")

    first = True
    with open("${uk_variants}", "r") as COG_variants_in, \
         open("${gisaid_variants}", "r") as GISAID_variants_in, \
         open("cog_gisaid_variants.csv", "w") as variants_out:
        for line in GISAID_variants_in:
            if first:
                variants_out.write(line)
                first = False
                continue

            sample = line.strip().split(",")[0]
            if sample in GISAID_fasta:
                variants_out.write(line)

        first = True
        for line in COG_variants_in:
            if first:
                first = False
                continue

            sample = line.strip().split(",")[0]
            if sample in COG_fasta:
                variants_out.write(line)
    """
}


process uk_geography {
    /**
    * Cleans up geography
    * @input uk_fasta, uk_metadata
    * @output geography_metadata
    * @params geography_utils
    */

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


process add_geography_to_metadata {
    /**
    * Adds UK geography to combined metadata
    * @input combined_metadata, geography_metadata
    * @output metadata
    */

    publishDir "${publish_dev}/COG_GISAID", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_master.csv"}

    input:
    path combined_metadata
    path geography_metadata

    output:
    path "cog_gisaid_geography.csv", emit: metadata

    script:
    """
    fastafunk add_columns \
          --in-metadata ${combined_metadata} \
          --in-data ${geography_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns adm1 adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
          --out-metadata "cog_gisaid_geography.csv"
    """
}


process split_recipes {
    output:
    path "*.json"

    script:
    """
    #!/usr/bin/env python3
    import json
    i = 0

    with open("${recipes}", 'r') as f:
        recipes = json.load(f)

        for d in recipes:
            for entry in recipes[d]:
                new_recipes = {d:[entry]}
                with open("%i.json" %i, 'w') as handle:
                    json.dump(new_recipes,handle)
                i += 1
    """
}


process publish_recipes {
    /**
    * Publishes subsets of combined FASTA and METADATA for COG-UK and GISAID
    * @input uk_unaligned_fasta, uk_aligned_fasta, uk_trimmed_fasta, combined_fasta,
    * uk_metadata, combined_metadata, uk_variants, combined_variants
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    input:
    path uk_unaligned_fasta
    path uk_aligned_fasta
    path uk_trimmed_fasta
    path combined_fasta
    path uk_metadata
    path combined_metadata
    path uk_variants
    path combined_variants
    path recipe

    output:
    path "*/cog_*.*", emit: csv

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --unaligned_fasta ${uk_unaligned_fasta} \
      --aligned_fasta ${uk_aligned_fasta} \
      --trimmed_fasta ${uk_trimmed_fasta} \
      --cog_global_fasta ${combined_fasta} \
      --cog_metadata ${uk_metadata} \
      --cog_global_metadata ${combined_metadata} \
      --cog_variants ${uk_variants} \
      --cog_global_variants ${combined_variants} \
      --recipes ${recipe} \
      --date ${params.date}
    """
}


process announce_to_webhook {
    input:
    path published_files

    script:
    if (params.webhook)
        """
        echo '{{"text":"' > announce.json
        echo "*Datapipe Complete*\\n" >> announce.json
        echo "> Dev outputs in : ${publish_dev}\\n" >> announce.json
        ls -R ${publish_dev} >> announce.json
        echo "> Publishable outputs in : ${publish_dir}\\n" >> announce.json
        ls -R ${publish_dir} >> announce.json
        echo '"}}' >> announce.json
        echo 'webhook ${params.webhook}'

        curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
        """
    else
        """
        touch "announce.json"
        """
}


geography_utils = file(params.uk_geography)
recipes = file(params.publish_recipes)

workflow publish_all {
    take:
        uk_unaligned_fasta
        uk_aligned_fasta
        uk_fasta
        uk_metadata
        uk_variants
        gisaid_fasta
        gisaid_metadata
        gisaid_variants
    main:
        combine_cog_gisaid(uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata)
        combine_variants(uk_fasta, uk_variants, gisaid_fasta, gisaid_variants)
        uk_geography(uk_fasta, uk_metadata)
        add_geography_to_metadata(combine_cog_gisaid.out.metadata,uk_geography.out.geography)
        split_recipes()
        recipe_ch = split_recipes.out.flatten()
        publish_recipes(uk_unaligned_fasta,uk_aligned_fasta,uk_fasta,combine_cog_gisaid.out.fasta, \
                        uk_metadata,add_geography_to_metadata.out.metadata,uk_variants,combine_variants.out, \
                        recipe_ch)
        outputs_ch = publish_recipes.out.csv.collect()
        announce_to_webhook(outputs_ch)
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    uk_variants = file(params.uk_variants)

    gisaid_fasta = file(params.gisaid_fasta)
    gisaid_metadata = file(params.gisaid_metadata)
    gisaid_variants = file(params.gisaid_variants)

    publish_all(uk_fasta,
                   uk_fasta,
                   uk_fasta,
                   uk_metadata,
                   uk_variants,
                   gisaid_fasta,
                   gisaid_metadata,
                   gisaid_variants)
}
