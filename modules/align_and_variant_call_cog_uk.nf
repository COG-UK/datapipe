#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process uk_minimap2_to_reference {
    /**
    * Minimaps samples to reference
    * @input uk_fasta
    * @output uk_sam
    * @params reference_fasta
    */

    cpus 8

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

    cpus 8

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

    cpus 8

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

    publishDir "${publish_dev}/COG", pattern: "*.csv", mode: 'copy', saveAs: {"cog_variants.csv"}

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
      --dels ${dels} \
      --index-column query
    """
}

workflow align_and_variant_call_cog_uk {
    take:
        uk_fasta
    main:
        uk_minimap2_to_reference(uk_fasta)
        uk_get_variants(uk_minimap2_to_reference.out)
        uk_get_indels(uk_minimap2_to_reference.out)
        uk_alignment(uk_minimap2_to_reference.out)
        uk_mask_alignment(uk_alignment.out)
        uk_get_snps(uk_mask_alignment.out)
        uk_type_AAs_and_dels(uk_mask_alignment.out, uk_get_variants.out)
    emit:
        variants = uk_type_AAs_and_dels.out
        fasta = uk_mask_alignment.out
}


mask_file = file(params.mask_file)
aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)


workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)

    align_and_variant_call_cog_uk(uk_fasta)
}