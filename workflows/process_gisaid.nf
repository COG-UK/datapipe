#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { preprocess_gisaid } from '../modules/preprocess_gisaid.nf'
include { pangolin } from '../modules/pangolin.nf'
include { check_for_pangolin_update } from '../modules/pangolin.nf'
include { deduplicate_gisaid } from '../modules/deduplicate.nf'
include { align_and_variant_call } from '../modules/align_and_variant_call.nf'
include { filter_and_trim_gisaid } from '../modules/filter_and_trim.nf'
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
      align_and_variant_call(deduplicate_gisaid.out.fasta, "gisaid")
      filter_and_trim_gisaid(align_and_variant_call.out.fasta, deduplicate_gisaid.out.metadata)
      publish_gisaid(filter_and_trim_gisaid.out.fasta, filter_and_trim_gisaid.out.metadata, align_and_variant_call.out.variants)
    emit:
      fasta = publish_gisaid.out.fasta
      metadata = publish_gisaid.out.metadata
      variants = publish_gisaid.out.variants
}

workflow {
    ch_gisaid_json = Channel.fromPath(params.gisaid_json)

    check_for_pangolin_update()
    process_gisaid(ch_gisaid_json, check_for_pangolin_update.out)
    announce_to_webhook(process_gisaid.out.published)

}
