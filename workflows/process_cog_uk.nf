#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { preprocess_cog_uk } from '../modules/preprocess_cog_uk.nf'
include { pangolin_cog_uk } from '../modules/pangolin_cog_uk.nf'
include { deduplicate_cog_uk } from '../modules/deduplicate_cog_uk.nf'
include { align_and_variant_call_cog_uk } from '../modules/align_and_variant_call_cog_uk.nf'
include { filter_and_trim_cog_uk } from '../modules/filter_and_trim_cog_uk.nf'
include { publish_all } from '../modules/publish_all.nf'

workflow process_cog_uk {
    take:
      uk_fasta
      uk_metadata
      uk_accessions
    main:
      preprocess_cog_uk(uk_fasta, uk_metadata, uk_accessions)
      pangolin_cog_uk(preprocess_cog_uk.out.fasta, preprocess_cog_uk.out.metadata)
      deduplicate_by_cogid_cog_uk(preprocess_cog_uk.out.fasta, pangolin_cog_uk.out.metadata)
      align_and_variant_call_cog_uk(deduplicate_cog_uk.out.fasta)
      deduplicate_by_biosample_cog_uk(align_and_variant_call_cog_uk.out.fasta,deduplicate_by_cogid_cog_uk.out.metadata)
      filter_and_trim_cog_uk(deduplicate_by_biosample_cog_uk.out.fasta, deduplicate_by_biosample_cog_uk.out.metadata)
    emit:
      unaligned_fasta = deduplicate_by_cogid_cog_uk.out.fasta
      aligned_fasta = align_and_variant_call_cog_uk.out.fasta
      trimmed_fasta = filter_and_trim_cog_uk.out.fasta
      metadata = filter_and_trim_cog_uk.out.metadata
      variants = align_and_variant_call_cog_uk.out.variants
}

workflow {
    ch_uk_fasta = Channel.fromPath(params.uk_fasta)
    ch_uk_metadata = Channel.fromPath(params.uk_metadata)
    ch_uk_accessions = Channel.fromPath(params.uk_accessions)

    process_cog_uk(ch_uk_fasta,
                   ch_uk_metadata,
                   ch_uk_accessions)

    ch_gisaid_fasta = Channel.fromPath(params.gisaid_fasta)
    ch_gisaid_metadata = Channel.fromPath(params.gisaid_metadata)
    ch_gisaid_variants = Channel.fromPath(params.gisaid_variants)

    publish_all(process_cog_uk.out.unaligned_fasta,
                process_cog_uk.out.aligned_fasta,
                process_cog_uk.out.trimmed_fasta,
                process_cog_uk.out.metadata,
                process_cog_uk.out.variants,
                ch_gisaid_fasta,
                ch_gisaid_metadata,
                ch_gisaid_variants)
}
