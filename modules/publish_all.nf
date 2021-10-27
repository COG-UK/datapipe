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
          --filter-column fasta_header covv_accession_id central_sample_id biosample_source_id secondary_identifier root_sample_id source_id \
                          sequence_name sample_date safe_sample_date epi_week epi_day collection_date received_date published_date \
                          country adm1 adm1_UK adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
                          is_uk is_cog_uk \
                          submission_org_code submission_user collection_pillar is_pillar_2 is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineages_version lineage_conflict lineage_ambiguity_score scorpio_call scorpio_support scorpio_conflict \
                          usher_lineage usher_lineages_version \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
                          unmapped_genome_completeness duplicate why_excluded nucleotide_mutations \
                          uk_lineage microreact_lineage del_lineage del_introduction phylotype \
          --where-column epi_week=edin_epi_week epi_day=edin_epi_day country=adm0 outer_postcode=adm2_private lineage_support=probability lineages_version=pangoLEARN_version adm1_UK=adm1 \
          --out-fasta "intermediate_cog.fa" \
          --out-metadata "intermediate_cog.csv" \
          --restrict --low-memory

        fastafunk fetch \
          --in-fasta ${gisaid_fasta} \
          --in-metadata ${gisaid_metadata} \
          --index-column sequence_name \
          --filter-column fasta_header covv_accession_id central_sample_id biosample_source_id secondary_identifier root_sample_id source_id \
                          sequence_name sample_date safe_sample_date epi_week epi_day collection_date received_date published_date \
                          country adm1 adm1_UK adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location \
                          is_uk is_cog_uk \
                          submission_org_code submission_user collection_pillar is_pillar_2 is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineages_version lineage_conflict lineage_ambiguity_score scorpio_call scorpio_support scorpio_conflict \
                          usher_lineage usher_lineages_version \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
                          unmapped_genome_completeness duplicate why_excluded nucleotide_mutations \
                          uk_lineage microreact_lineage del_lineage del_introduction phylotype \
          --where-column adm1=edin_admin_1 travel_history=edin_travel \
          --out-fasta "intermediate_gisaid.fa" \
          --out-metadata "intermediate_gisaid.csv" \
          --restrict --low-memory

        cat intermediate_cog.fa intermediate_gisaid.fa > cog_gisaid.fa
        cat intermediate_cog.csv > cog_gisaid.csv
        tail -n+2 intermediate_gisaid.csv >> cog_gisaid.csv

        head -n1 intermediate_cog.csv > head_cog.txt
        head -n1 intermediate_gisaid.csv > head_gisaid.txt
        cmp --silent head_cog.txt head_gisaid.txt || exit 1
    """
}


process combine_mutations {
    /**
    * Combines FASTA and mutation metadata for COG-UK and GISAID
    * @input uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata
    * @output cog_gisaid_fasta, cog_gisaid_metadata
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_mutations.csv"}

    input:
    path uk_mutations
    path gisaid_mutations

    output:
    path "cog_gisaid_mutations.csv"

    script:
    """
    fastafunk merge \
      --in-metadata ${uk_mutations} ${gisaid_mutations} \
      --out-metadata "cog_gisaid_mutations.csv" \
      --index-column "sequence_name"
    """
}

process combine_constellations {
    /**
    * Combines FASTA and constellation metadata for COG-UK and GISAID
    * @input uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata
    * @output cog_gisaid_fasta, cog_gisaid_metadata
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_constellations.csv"}

    input:
    path uk_constellations
    path gisaid_constellations

    output:
    path "cog_gisaid_constellations.csv"

    script:
    """
    fastafunk merge \
      --in-metadata ${uk_constellations} ${gisaid_constellations} \
      --out-metadata "cog_gisaid_constellations.csv" \
      --index-column "sequence_name"
    """
}

process combine_updown {
    /**
    * Combines updown metadata for COG-UK and GISAID
    * @input uk_updown gisaid_updown
    * @output cog_gisaid_updown
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_updown.csv"}

    input:
    path uk_updown
    path gisaid_updown

    output:
    path "cog_gisaid_updown.csv"

    script:
    """
    cp ${uk_updown} tmp.csv
    tail -n+1 ${gisaid_updown} >> tmp.csv
    grep -v ",,,," tmp.csv > "cog_gisaid_updown.csv"
    """
}


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
      --filter-column central_sample_id sequence_name sample_date epi_week \
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
      --epiweek-col epi_week \
      --outdir geography

    rm -rf geography_tmp
    """
}


process add_uk_geography_to_metadata {
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
          --force-overwrite \
          --new-columns adm1 adm1_raw adm2 outer_postcode adm2_raw adm2_source NUTS1 region latitude longitude location safe_location utla utla_code suggested_adm2_grouping \
          --out-metadata "cog_gisaid_geography.csv"
    """
}


process gisaid_geography {
    /**
    * Cleans up geography
    * @input fasta, metadata
    * @output geography_metadata
    * @params geography_utils
    */

    memory { 1.GB * task.attempt + fasta.size() * 1.B }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries = 1

    publishDir "${publish_dev}/", pattern: "geography/*.csv", mode: 'copy'
    publishDir "${publish_dev}/", pattern: "geography/*.txt", mode: 'copy'

    input:
    path fasta
    path metadata

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
    * @input combined_metadata, geography_metadata
    * @output metadata
    */

    publishDir "${publish_dev}/gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"gisaid_master.csv"}, overwrite: true

    input:
    path combined_metadata
    path geography_metadata

    output:
    path "gisaid_geography.csv", emit: metadata

    script:
    """
    fastafunk add_columns \
          --in-metadata ${combined_metadata} \
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
    * uk_metadata, combined_metadata, uk_mutations, combined_mutations
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false
    publishDir "${publish_dir}/", pattern: "README", mode: 'copy', overwrite: false

    memory {8.GB * task.attempt}
    errorStrategy = { 'retry' }
    maxRetries 3

    input:
    tuple path(uk_unaligned_fasta),path(uk_aligned_fasta),path(uk_trimmed_fasta),path(combined_fasta),path(uk_metadata),path(combined_metadata),path(combined_mutations),path(combined_constellations),path(combined_updown),path(recipe)

    output:
    path "${recipe.baseName}.done.txt", emit: flag
    path "README", emit: readme
    path "public/cog_${params.date}_all.fa", optional: true, emit: fasta
    path "public/cog_${params.date}_metadata.csv", optional: true, emit: metadata
    path "public/cog_${params.date}_alignment.fa", optional: true, emit: alignment
    path "public/cog_${params.date}_unmasked_alignment.fa", optional: true, emit: unmasked_alignment
    path "*/cog_*.*", emit: all

    script:
    """
    cp $project_dir/../resources/publish_readme.txt README

    $project_dir/../bin/publish_from_config.py \
      --unaligned_fasta ${uk_unaligned_fasta} \
      --aligned_fasta ${uk_aligned_fasta} \
      --trimmed_fasta ${uk_trimmed_fasta} \
      --cog_global_fasta ${combined_fasta} \
      --cog_metadata ${uk_metadata} \
      --cog_global_metadata ${combined_metadata} \
      --mutations ${combined_mutations} \
      --constellations ${combined_constellations} \
      --updown ${combined_updown} \
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
    publishDir "${publish_dev}/", pattern: "s3dir", mode: 'copy'

    input:
    path fasta
    path metadata
    path alignment
    path unmasked_alignment

    output:
    path s3dir


    script:
    """
    mkdir -p s3dir
    cp ${fasta} s3dir/cog_all.fasta
    cp ${metadata} s3dir/cog_metadata.csv
    cp ${alignment} s3dir/cog_alignment.fasta
    cp ${unmasked_alignment} s3dir/cog_unmasked_alignment.fasta
    """
}


process publish_gisaid_recipes {
    /**
    * Publishes subsets of combined FASTA and METADATA for COG-UK and GISAID
    * @input gisaid_unaligned_fasta, gisaid_aligned_fasta, gisaid_trimmed_fasta, combined_fasta,
    * gisaid_metadata, combined_metadata, gisaid_mutations, combined_mutations
    * @params publish_recipes.json
    * @output many
    */

    publishDir "${publish_dir}/", pattern: "*/*.*", mode: 'copy', overwrite: false

    memory {8.GB * task.attempt}
    errorStrategy = { 'retry' }
    maxRetries 3

    input:
    tuple path(gisaid_fasta),path(gisaid_metadata),path(gisaid_mutations),path(gisaid_constellations),path(gisaid_updown),path(recipe)

    output:
    path "*/gisaid_*.*", emit: all
    path "*/gisaid_*_global_alignment.fa", optional: true, emit: fasta
    path "*/gisaid_*_global_metadata.csv", optional: true, emit: metadata
    path "*/gisaid_*_global_mutations.csv", optional: true, emit: mutations
    path "*/gisaid_*_global_constellations.csv", optional: true, emit: constellations
    path "*/gisaid_*_global_updown.csv", optional: true, emit: updown

    script:
    """
    $project_dir/../bin/publish_from_config.py \
      --recipes ${recipe} \
      --date ${params.date} \
      --gisaid_fasta ${gisaid_fasta} \
      --gisaid_metadata ${gisaid_metadata} \
      --mutations ${gisaid_mutations} \
      --constellations ${gisaid_constellations} \
      --updown ${gisaid_updown}
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
        uk_mutations
        uk_constellations
        uk_updown
        gisaid_fasta
        gisaid_metadata
        gisaid_mutations
        gisaid_constellations
        gisaid_updown
    main:
        combine_cog_gisaid(uk_fasta, uk_metadata, gisaid_fasta, gisaid_metadata)
        combine_mutations(uk_mutations, gisaid_mutations)
        combine_constellations(uk_constellations, gisaid_constellations)
        combine_updown(uk_updown, gisaid_updown)
        uk_geography(uk_fasta, uk_metadata)
        add_uk_geography_to_metadata(combine_cog_gisaid.out.metadata,uk_geography.out.geography)
        make_delta_by_utla_summary(add_uk_geography_to_metadata.out.metadata)
        split_recipes(cog_global_recipes)
        recipe_ch = split_recipes.out.flatten()
        uk_unaligned_fasta.combine(uk_aligned_fasta)
                          .combine(uk_fasta)
                          .combine(combine_cog_gisaid.out.fasta)
                          .combine(uk_metadata)
                          .combine(add_uk_geography_to_metadata.out.metadata)
                          .combine(combine_mutations.out)
                          .combine(combine_constellations.out)
                          .combine(combine_updown.out)
                          .combine(recipe_ch)
                          .set{ publish_input_ch }
        publish_cog_global_recipes(publish_input_ch)
        outputs_ch = publish_cog_global_recipes.out.flag.collect()
        announce_to_webhook(outputs_ch, "${params.whoami}")
        if ( params.s3 )
        {
            publish_s3(publish_cog_global_recipes.out.fasta, publish_cog_global_recipes.out.metadata, publish_cog_global_recipes.out.alignment, publish_cog_global_recipes.out.unmasked_alignment)
        }
}


workflow publish_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
        gisaid_mutations
        gisaid_constellations
        gisaid_updown
    main:
        if ( params.geography ){
            gisaid_geography(gisaid_fasta, gisaid_metadata)
            add_gisaid_geography_to_metadata(gisaid_metadata,gisaid_geography.out.geography)
            add_gisaid_geography_to_metadata.out.metadata.set{ new_gisaid_metadata }
        } else {
            new_gisaid_metadata = gisaid_metadata
        }
        split_recipes(gisaid_recipes)
        recipe_ch = split_recipes.out.flatten()
        gisaid_fasta.combine(new_gisaid_metadata)
                    .combine(gisaid_mutations)
                    .combine(gisaid_constellations)
                    .combine(gisaid_updown)
                    .combine(recipe_ch)
                    .set{ publish_input_ch }
        publish_gisaid_recipes(publish_input_ch)
        outputs_ch = publish_gisaid_recipes.out.all.collect()
    emit:
        fasta = publish_gisaid_recipes.out.fasta
        metadata = publish_gisaid_recipes.out.metadata
        mutations = publish_gisaid_recipes.out.mutations
        constellations = publish_gisaid_recipes.out.constellations
        updown = publish_gisaid_recipes.out.updown
        published = outputs_ch
}


workflow {
    uk_unaligned_fasta = Channel.fromPath(params.uk_unaligned_fasta)
    uk_aligned_fasta = Channel.fromPath(params.uk_aligned_fasta)
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    uk_mutations = Channel.fromPath(params.uk_mutations)
    uk_constellations = Channel.fromPath(params.uk_constellations)
    uk_updown = Channel.fromPath(params.uk_updown)

    gisaid_fasta = Channel.fromPath(params.gisaid_fasta)
    gisaid_metadata = Channel.fromPath(params.gisaid_metadata)
    gisaid_mutations = Channel.fromPath(params.gisaid_mutations)
    gisaid_constellations = Channel.fromPath(params.gisaid_constellations)
    gisaid_updown = Channel.fromPath(params.gisaid_updown)

    publish_all(uk_unaligned_fasta,
                   uk_aligned_fasta,
                   uk_fasta,
                   uk_metadata,
                   uk_mutations,
                   uk_constellations,
                   uk_updown,
                   gisaid_fasta,
                   gisaid_metadata,
                   gisaid_mutations,
                   gisaid_constellations,
                   gisaid_updown)
}
