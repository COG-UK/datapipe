#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { preprocess_gisaid } from '../modules/preprocess_gisaid.nf'
include { pangolin } from '../modules/pangolin.nf'
include { deduplicate_gisaid } from '../modules/deduplicate.nf'
include { align_and_variant_call } from '../modules/align_and_variant_call.nf'
include { filter_and_trim_gisaid } from '../modules/filter_and_trim.nf'
include { publish_gisaid } from '../modules/publish_all.nf'

workflow process_gisaid {
    take:
      gisaid_json
    main:
      preprocess_gisaid(gisaid_json)
      pangolin(preprocess_gisaid.out.fasta, preprocess_gisaid.out.metadata)
      deduplicate_gisaid(preprocess_gisaid.out.fasta, pangolin.out.metadata)
      align_and_variant_call(deduplicate_gisaid.out.fasta)
      filter_and_trim_gisaid(align_and_variant_call.out.fasta, deduplicate_gisaid.out.metadata)
      publish_gisaid(filter_and_trim_gisaid.out.fasta, filter_and_trim_gisaid.out.metadata, align_and_variant_call.out.variants)
    emit:
      fasta = filter_and_trim_gisaid.out.fasta
      metadata = filter_and_trim_gisaid.out.metadata
      variants = align_and_variant_call.out.variants
}

workflow {
    ch_gisaid_json = Channel.fromPath(params.gisaid_json)

    process_gisaid(ch_gisaid_json)

}
