#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include { preprocess_cog_uk } from '../modules/preprocess_cog_uk.nf'
include { pangolin_cog_uk } from '../modules/pangolin_cog_uk.nf'
include { deduplicate_cog_uk } from '../modules/deduplicate_cog_uk.nf'
include { align_and_variant_call_cog_uk } from '../modules/align_and_variant_call_cog_uk.nf'
include { filter_and_trim_cog_uk } from '../modules/filter_and_trim_cog_uk.nf'
include { publish_cog_uk } from '../modules/publish_cog_uk.nf'

workflow process_cog_uk {
    take:
      uk_fasta
      uk_metadata
    main:
      preprocess_cog_uk(uk_fasta, uk_metadata)
      pangolin_cog_uk(preprocess_cog_uk.out.fasta, preprocess_cog_uk.out.metadata)
      deduplicate_cog_uk(preprocess_cog_uk.out.fasta, pangolin_cog_uk.out.metadata)
      align_and_variant_call_cog_uk(deduplicate_cog_uk.out.fasta, deduplicate_cog_uk.out.metadata)
      filter_and_trim_cog_uk(align_and_variant_call_cog_uk.out.fasta, align_and_variant_call_cog_uk.out.metadata)
    emit:
}