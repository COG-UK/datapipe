#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)


process uk_minimap2_to_reference {
    /**
    * Minimaps samples to reference
    * @input uk_fasta
    * @output uk_sam
    * @params reference_fasta
    */

    cpus 1

    input:
    path uk_fasta

    output:
    path "uk_alignment.sam"

    script:
    """
    minimap2 -t ${task.cpus} -a -x asm5 ${reference_fasta} ${uk_fasta} > uk_alignment.sam
    """
}

process uk_get_variants {
    /**
    * Creates CSV of variants found in each genome
    * @input uk_sam
    * @output uk_variants
    * @parms reference_fasta, reference_genbank
    */

    cpus 1

    input:
    path uk_sam

    output:
    path "uk.variants.csv"

    script:
    """
    gofasta sam variants -t ${task.cpus} \
      --samfile ${uk_sam} \
      --reference ${reference_fasta} \
      --genbank ${reference_genbank} \
      --outfile uk.variants.csv
    """
}

process uk_get_indels {
    /**
    * Creates CSV of indels found in each genome
    * @input uk_sam
    * @output uk_insertions, uk_deletions
    */

    publishDir "${publish_dir}/metadata/", pattern: "*.csv", mode: 'copy'

    input:
    path uk_sam

    output:
    path "uk.insertions.csv", emit: uk_insertions
    path "uk.deletions.csv", emit: uk_deletions

    script:
    """
    gofasta sam indels \
      -s ${uk_sam} \
      --threshold 2 \
      --insertions-out "uk.insertions.csv" \
      --deletions-out "uk.deletions.csv"
    """
}

process uk_alignment {
    /**
    * Get reference-based alignment
    * @input uk_sam
    * @output uk_alignment
    * @params reference_fasta
    */

    cpus 1

    input:
    path uk_sam

    output:
    path "uk_alignment.fasta"

    script:
    """
    gofasta sam toMultiAlign -t ${task.cpus} \
      --samfile ${uk_sam} \
      --reference ${reference_fasta} \
      -o uk_alignment.fasta
    """
}


process uk_mask_alignment {
    /**
    * Applies a mask to aligned FASTA
    * @input uk_alignment
    * @output uk_alignment_updated
    * @params mask_file
    */

    input:
    path uk_alignment

    output:
    path "${uk_alignment.baseName}.masked.fa"

    script:
    """
    $project_dir/../bin/add_mask.py \
      --in-alignment ${uk_alignment} \
      --out-alignment "${uk_alignment.baseName}.masked.fa" \
      --mask ${mask_file} \
    """
}


process uk_get_snps {
    /**
    * Call SNPs in each genome
    * @input uk_alignment
    * @output uk_snps
    * @params reference_fasta
    */

    input:
    path uk_alignment

    output:
    path "uk.snps.csv"

    script:
    """
    gofasta snps -r ${reference_fasta} -q ${uk_alignment} -o uk.snps.csv
    """
}

process uk_type_AAs_and_dels {
    /**
    * Adds a column to metadata table for specific dels and aas looked for
    * @input uk_alignment, uk_metadata
    * @output uk_metadata_updated
    * @params reference_fasta, del_file, aa_file
    */

    input:
    path uk_alignment
    path uk_metadata

    output:
    path "${uk_metadata.baseName}.typed.csv"

    script:
    """
    $project_dir/../bin/type_aas_and_dels.py \
      --in-fasta ${uk_alignment} \
      --in-metadata ${uk_metadata} \
      --out-metadata ${uk_metadata.baseName}.typed.csv \
      --reference-fasta ${reference_fasta} \
      --aas ${aas} \
      --dels ${dels}
    """
}

process publish_full_aligned_cog_data {
    /**
    * Publish full alignment and relevant metadata
    * @input uk_alignment, uk_metadata
    * @params date
    */

    publishDir "${publish_dir}/alignments/", pattern: "*.fa", mode: 'copy', saveAs: {"cog_${params.date}_all_alignment.fasta"}
    publishDir "${publish_dir}/alignments/", pattern: "*.csv", mode: 'copy', saveAs: {"cog_${params.date}_all_metadata.fasta"}


    input:
    path uk_alignment
    path uk_metadata

    output:
    path "${uk_alignment.baseName}.matched.fa"
    path "${uk_metadata.baseName}.matched.csv"

    script:
    """
        fastafunk fetch \
          --in-fasta ${uk_alignment} \
          --in-metadata ${uk_metadata} \
          --index-column sequence_name \
          --filter-column \
                country adm1 adm2 outer_postcode biosample_source_id \
                central_sample_id collected_by collection_date end_time \
                flowcell_id flowcell_type instrument_make instrument_model is_surveillance \
                layout_insert_length layout_read_length \
                library_adaptor_barcode library_layout_config library_name library_primers library_protocol \
                library_selection library_seq_kit library_seq_protocol library_source library_strategy \
                meta.artic.primers meta.artic.protocol meta.epi.cluster meta.investigation.cluster \
                meta.investigation.name meta.investigation.site metric.ct.1.ct_value metric.ct.1.test_kit \
                metric.ct.1.test_platform metric.ct.1.test_target metric.ct.2.ct_value metric.ct.2.test_kit         \
                metric.ct.2.test_platform metric.ct.2.test_target metric.ct.max_ct \
                metric.ct.min_ct metric.ct.num_tests \
                published_as received_date root_sample_id run_group run_name \
                sample_type_collected sample_type_received secondary_accession secondary_identifier \
                sequencing_org sequencing_org_code sequencing_submission_date sequencing_uuid \
                source_age source_sex start_time \
                submission_org submission_org_code submission_user swab_site \
                header sequence_name length missing gaps cov_id sample_date subsample_omit epi_week \
                d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3 \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private \
          --out-fasta "${uk_alignment.baseName}.matched.fa" \
          --out-metadata "${uk_metadata.baseName}.matched.csv" \
          --restrict
        """
}

workflow align_and_variant_call_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        uk_minimap2_to_reference(uk_fasta)
        uk_get_variants(uk_minimap2_to_reference.out)
        uk_get_indels(uk_minimap2_to_reference.out)
        uk_alignment(uk_minimap2_to_reference.out)
        uk_mask_alignment(uk_alignment.out)
        uk_get_snps(uk_mask_alignment.out)
        uk_type_AAs_and_dels(uk_mask_alignment.out, uk_metadata)
        publish_full_aligned_cog_data(uk_mask_alignment.out, uk_type_AAs_and_dels.out)
    emit:
        variants = uk_get_variants.out
        fasta = uk_alignment.out
        metadata = uk_type_AAs_and_dels.out
}


mask_file = file(params.mask_file)
aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)


workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)

    align_and_variant_call_cog_uk(uk_fasta,
                                  uk_metadata)
}