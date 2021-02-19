#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import modules
include { check_for_pangolin_update } from '../modules/pangolin.nf'
include { process_gisaid } from './process_gisaid.nf'
include { process_cog_uk } from './process_cog_uk.nf'
include { publish_cog_global } from './process_cog_uk.nf'


workflow {
    ch_uk_fasta = Channel.fromPath(params.uk_fasta)
    ch_uk_metadata = Channel.fromPath(params.uk_metadata)
    ch_uk_accessions = Channel.fromPath(params.uk_accessions)
    ch_gisaid_json = Channel.fromPath(params.gisaid_json)

    check_for_pangolin_update()
    process_gisaid(ch_gisaid_json, check_for_pangolin_update.out)
    process_cog_uk(ch_uk_fasta,
                   ch_uk_metadata,
                   ch_uk_accessions,
                   check_for_pangolin_update.out)

    publish_cog_global(process_cog_uk.out.unaligned_fasta,
                        process_cog_uk.out.aligned_fasta,
                        process_cog_uk.out.trimmed_fasta,
                        process_cog_uk.out.metadata,
                        process_cog_uk.out.variants,
                        process_gisaid.out.fasta,
                        process_gisaid.out.metadata,
                        process_gisaid.out.variants)
}
