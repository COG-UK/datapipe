#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { start } from '../modules/start.nf'
include { preprocess_cog_uk } from '../modules/preprocess_cog_uk.nf'
include { pangolin } from '../modules/pangolin.nf'
include { check_for_pangolin_update } from '../modules/pangolin.nf'
include { deduplicate_cog_uk } from '../modules/deduplicate.nf'
include { align_and_variant_call } from '../modules/align_and_variant_call.nf'
include { filter_and_trim_cog_uk } from '../modules/filter_and_trim.nf'
include { publish_cog_global } from '../modules/publish_all.nf'

workflow process_cog_uk {
    take:
      uk_fasta
      uk_metadata
      uk_accessions
      pangolin_updated
    main:
      preprocess_cog_uk(uk_fasta, uk_metadata, uk_accessions)
      pangolin(preprocess_cog_uk.out.fasta, preprocess_cog_uk.out.metadata, pangolin_updated)
      deduplicate_cog_uk(preprocess_cog_uk.out.fasta, pangolin.out.metadata)
      align_and_variant_call(deduplicate_cog_uk.out.fasta, deduplicate_cog_uk.out.metadata, "cog")
      filter_and_trim_cog_uk(align_and_variant_call.out.fasta, align_and_variant_call.out.metadata)
    emit:
      unaligned_fasta = deduplicate_cog_uk.out.fasta
      aligned_fasta = align_and_variant_call.out.fasta
      trimmed_fasta = filter_and_trim_cog_uk.out.fasta
      metadata = filter_and_trim_cog_uk.out.metadata
      mutations = align_and_variant_call.out.mutations
      constellations = align_and_variant_call.out.constellations
      updown = align_and_variant_call.out.updown
}

workflow {
    start()

    ch_uk_fasta = Channel.fromPath(params.uk_fasta)
    ch_uk_metadata = Channel.fromPath(params.uk_metadata)
    ch_uk_accessions = Channel.fromPath(params.uk_accessions)

    check_for_pangolin_update()
    process_cog_uk(ch_uk_fasta,
                   ch_uk_metadata,
                   ch_uk_accessions,
                   check_for_pangolin_update.out)

    ch_gisaid_fasta = Channel.fromPath(params.gisaid_fasta)
    ch_gisaid_metadata = Channel.fromPath(params.gisaid_metadata)
    ch_gisaid_mutations = Channel.fromPath(params.gisaid_mutations)
    ch_gisaid_constellations = Channel.fromPath(params.gisaid_constellations)

    publish_cog_global(process_cog_uk.out.unaligned_fasta,
                        process_cog_uk.out.aligned_fasta,
                        process_cog_uk.out.trimmed_fasta,
                        process_cog_uk.out.metadata,
                        process_cog_uk.out.mutations,
                        process_cog_uk.out.constellations,
                        ch_gisaid_fasta,
                        ch_gisaid_metadata,
                        ch_gisaid_mutations,
                        ch_gisaid_constellations)
}
