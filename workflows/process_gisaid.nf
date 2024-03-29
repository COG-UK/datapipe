#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { start } from '../modules/start.nf'
include { preprocess_gisaid } from '../modules/preprocess_gisaid.nf'
include { pangolin } from '../modules/pangolin.nf'
include { check_for_pangolin_update } from '../modules/pangolin.nf'
include { deduplicate_gisaid } from '../modules/deduplicate.nf'
include { align_and_variant_call } from '../modules/align_and_variant_call.nf'
include { filter_and_trim_gisaid } from '../modules/filter_and_trim.nf'
include { clean_geography_gisaid } from '../modules/clean_geography.nf'
include { publish_gisaid } from '../modules/publish_all.nf'
include { announce_to_webhook } from '../modules/publish_all.nf'

workflow process_gisaid {
    take:
      gisaid_json
      pangolin_updated
    main:
      preprocess_gisaid(gisaid_json)
      pangolin(preprocess_gisaid.out.fasta, preprocess_gisaid.out.metadata, pangolin_updated)
      deduplicate_gisaid(preprocess_gisaid.out.fasta, pangolin.out.metadata)
      align_and_variant_call(deduplicate_gisaid.out.fasta, deduplicate_gisaid.out.metadata, "gisaid")
      filter_and_trim_gisaid(align_and_variant_call.out.fasta, align_and_variant_call.out.metadata)
      clean_geography_gisaid(filter_and_trim_gisaid.out.fasta, filter_and_trim_gisaid.out.metadata)
      publish_gisaid(filter_and_trim_gisaid.out.fasta, clean_geography_gisaid.out.metadata, align_and_variant_call.out.mutations, align_and_variant_call.out.constellations, align_and_variant_call.out.updown)
    emit:
      fasta = publish_gisaid.out.fasta
      metadata = publish_gisaid.out.metadata
      mutations = publish_gisaid.out.mutations
      constellations = publish_gisaid.out.constellations
      updown = publish_gisaid.out.updown
}

workflow {
    start()
    
    ch_gisaid_json = Channel.fromPath(params.gisaid_json)

    check_for_pangolin_update()
    process_gisaid(ch_gisaid_json, check_for_pangolin_update.out)
    announce_to_webhook(process_gisaid.out.metadata, "Gisaid processing")

}
