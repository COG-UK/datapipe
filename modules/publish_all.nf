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

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.fa", mode: 'copy'

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
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
                          unmapped_genome_completeness duplicate why_excluded nucleotide_variants \
                          uk_lineage microreact_lineage del_lineage del_introduction phylotype \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private lineage_support=probability lineages_version=pangoLEARN_version adm1_UK=adm1\
          --out-fasta "intermediate_cog.fa" \
          --out-metadata "intermediate_cog.csv" \
          --restrict --low-memory

        fastafunk fetch \
          --in-fasta ${gisaid_fasta} \
          --in-metadata ${gisaid_metadata} \
          --index-column sequence_name \
          --filter-column fasta_header covv_accession_id central_sample_id biosample_source_id secondary_identifier root_sample_id \
                          pillar_2 \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
                          submission_org_code is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineage_support lineages_version \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
                          unmapped_genome_completeness duplicate why_excluded nucleotide_variants\
                          uk_lineage microreact_lineage del_lineage del_introduction phylotype \
          --where-column adm1=edin_admin_1 travel_history=edin_travel \
          --out-fasta "intermediate_gisaid.fa" \
          --out-metadata "intermediate_gisaid.csv" \
          --restrict --low-memory

        fastafunk merge \
          --in-fasta "intermediate_cog.fa" "intermediate_gisaid.fa" \
          --in-metadata "intermediate_cog.csv" "intermediate_gisaid.csv" \
          --out-fasta "cog_gisaid.fa" \
          --out-metadata "cog_gisaid.csv" \
          --index-column sequence_name \
          --low-memory
    """
}


process combine_variants {
    /**
    * Combines FASTA and METADATA for COG-UK and GISAID
    * @input uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata
    * @output cog_gisaid_fasta, cog_gisaid_metadata
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_variants.csv"}


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
        for line in COG_variants_in:
            if first:
                variants_out.write(line)
                first = False
                continue

            sample = line.strip().split(",")[0]
            if sample in COG_fasta:
                variants_out.write(line)

        first = True
        for line in GISAID_variants_in:
            if first:
                first = False
                continue

            sample = line.strip().split(",")[0]
            if sample in GISAID_fasta:
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

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_master.csv"}

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
    input:
    path recipes

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


process publish_cog_global_recipes {
    /**
    * Publishes subsets of combined FASTA and METADATA for COG-UK and GISAID
    * @input uk_unaligned_fasta, uk_aligned_fasta, uk_trimmed_fasta, combined_fasta,
    * uk_metadata, combined_metadata, uk_variants, combined_variants
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    input:
    tuple path(uk_unaligned_fasta),path(uk_aligned_fasta),path(uk_trimmed_fasta),path(combined_fasta),path(uk_metadata),path(combined_metadata),path(uk_variants),path(combined_variants),path(recipe)

    output:
    path "${recipe.baseName}.done.txt", emit: all
    path "public/cog_*_all.fa", optional: true, emit: fasta
    path "public/cog_*_metadata.csv", optional: true, emit: metadata
    path "public/cog_*_alignment.fa", optional: true, emit: alignment
    path "public/cog_*_unmasked_alignment.fa", optional: true, emit: unmasked_alignment

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
      touch "${recipe.baseName}.done.txt"
    """
}

process publish_s3 {
    /**
    * Publishes public files to s3
    * @input fasta, metadata, aligment, unmasked_alignment
    */

    input:
    path fasta
    path metadata
    path alignment
    path unmasked_alignment

    script:
    """
    mkdir -p s3dir
    cp ${fasta} s3dir/cog_all.fasta
    cp ${metadata} s3dir/cog_metadata.csv
    cp ${alignment} s3dir/cog_alignment.fasta
    cp ${unmasked_alignment} s3dir/cog_unmasked_alignment.fasta

    s3cmd sync s3dir/ s3://cog-uk/phylogenetics/{params.date}/ --acl-public
    s3cmd sync s3dir/ s3://cog-uk/phylogenetics/latest/ --acl-public
    """
}


process publish_gisaid_recipes {
    /**
    * Publishes subsets of combined FASTA and METADATA for COG-UK and GISAID
    * @input gisaid_unaligned_fasta, gisaid_aligned_fasta, gisaid_trimmed_fasta, combined_fasta,
    * gisaid_metadata, combined_metadata, gisaid_variants, combined_variants
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dev}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    input:
    tuple path(gisaid_fasta),path(gisaid_metadata),path(gisaid_variants),path(recipe)

    output:
    path "*/gisaid_*.*", emit: all
    path "*/gisaid_*_global_alignment.fa", optional: true, emit: fasta
    path "*/gisaid_*_global_metadata.csv", optional: true, emit: metadata
    path "*/gisaid_*_global_variants.csv", optional: true, emit: variants

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --recipes ${recipe} \
      --date ${params.date} \
      --gisaid_fasta ${gisaid_fasta} \
      --gisaid_metadata ${gisaid_metadata} \
      --gisaid_variants ${gisaid_variants}
    """
}


process announce_to_webhook {
    input:
    file published_files
    val name

    script:
    if (params.webhook)
        """
        echo '{"text":"' > announce.json
        echo "*${name} Complete*\\n" >> announce.json
        echo "> Dev outputs in : ${publish_dev}\\n" >> announce.json
        echo "> Publishable outputs in : ${publish_dir}\\n" >> announce.json
        echo '"}' >> announce.json
        echo 'webhook ${params.webhook}'

        curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
        """
    else
        """
        touch "announce.json"
        """
}


geography_utils = file(params.uk_geography)
cog_global_recipes = file(params.publish_cog_global_recipes)
gisaid_recipes = file(params.publish_gisaid_recipes)


workflow publish_cog_global {
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
        split_recipes(cog_global_recipes)
        recipe_ch = split_recipes.out.flatten()
        uk_unaligned_fasta.combine(uk_aligned_fasta)
                          .combine(uk_fasta)
                          .combine(combine_cog_gisaid.out.fasta)
                          .combine(uk_metadata)
                          .combine(add_geography_to_metadata.out.metadata)
                          .combine(uk_variants)
                          .combine(combine_variants.out)
                          .combine(recipe_ch)
                          .set{ publish_input_ch }
        publish_cog_global_recipes(publish_input_ch)
        outputs_ch = publish_cog_global_recipes.out.all.collect()
        announce_to_webhook(outputs_ch, "Datapipe")
        publish_s3(publish_cog_global_recipes.out.fasta, publish_cog_global_recipes.out.metadata, publish_cog_global_recipes.out.alignment, publish_cog_global_recipes.out.unmasked_alignment)
}


workflow publish_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
        gisaid_variants
    main:
        split_recipes(gisaid_recipes)
        recipe_ch = split_recipes.out.flatten()
        gisaid_fasta.combine(gisaid_metadata)
                    .combine(gisaid_variants)
                    .combine(recipe_ch)
                    .set{ publish_input_ch }
        publish_gisaid_recipes(publish_input_ch)
        outputs_ch = publish_gisaid_recipes.out.all.collect()
    emit:
        fasta = publish_gisaid_recipes.out.fasta
        metadata = publish_gisaid_recipes.out.metadata
        variants = publish_gisaid_recipes.out.variants
        published = outputs_ch
}


workflow {
    uk_unaligned_fasta = Channel.fromPath(params.uk_unaligned_fasta)
    uk_aligned_fasta = Channel.fromPath(params.uk_aligned_fasta)
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    uk_variants = Channel.fromPath(params.uk_variants)

    gisaid_fasta = Channel.fromPath(params.gisaid_fasta)
    gisaid_metadata = Channel.fromPath(params.gisaid_metadata)
    gisaid_variants = Channel.fromPath(params.gisaid_variants)

    publish_all(uk_unaligned_fasta,
                   uk_aligned_fasta,
                   uk_fasta,
                   uk_metadata,
                   uk_variants,
                   gisaid_fasta,
                   gisaid_metadata,
                   gisaid_variants)
}
