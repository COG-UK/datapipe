#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
includeConfig 'config/base.config'


process uk_minimap2_to_reference {
    /**
    * Minimaps samples to reference
    * @input uk_fasta
    * @output uk_sam
    * @params reference_fasta
    */

    cpus 1

    input:
    file uk_fasta

    output:
    path "uk_alignment.sam"

    script:
    """
    minimap2 -t ${task.cpus} -a -x asm5 ${params.reference_fasta} ${uk_fasta} > uk_alignment.sam
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
    file uk_sam

    output:
    path "uk.variants.csv"

    script:
    """
    gofasta sam variants -t ${task.cpus} \
      --samfile ${uk_sam} \
      --reference {params.reference_fasta} \
      --genbank {params.reference_genbank} \
      --outfile uk.variants.csv
    """
}

process uk_get_indels {
    /**
    * Creates CSV of indels found in each genome
    * @input uk_sam
    * @output uk_insertions, uk_deletions
    */

    input:
    file uk_sam

    output:
    path "uk.insertions.csv", emit: uk_insertions
    path "uk.deletions.csv", emit: uk_deletions

    script:
    """
    indels ${uk_sam} "uk.insertions.csv" "uk.deletions.csv"
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
    file uk_sam

    output:
    path "uk_alignment.fasta"

    script:
    """
    gofasta sam toMultiAlign -t ${task.cpus} \
      --samfile ${uk_sam} \
      --reference ${params.reference_fasta} \
      -o uk_alignment.fasta
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
    file uk_alignment

    output:
    path "uk.snps.csv"

    script:
    """
    snps ${params.reference_fasta} ${uk_alignment} > uk.snps.csv
    """
}

process type_AAs_and_dels {
    /**
    * Adds a column to metadata table for specific dels and aas looked for
    * @input uk_alignment, uk_metadata
    * @output uk_metadata_updated
    * @params reference_fasta, del_file, aa_file
    */

    input:
    file uk_alignment
    file uk_metadata

    output:
    path "${uk_metadata.baseName}.typed.csv"

    script:
    """
    $project_dir/../bin/type_aas_and_dels.py \
      --in-fasta ${uk_alignment} \
      --in-metadata ${uk_metadata} \
      --out-metadata ${uk_metadata.baseName}.typed.csv \
      --reference-fasta ${params.reference_fasta} \
      --aas ${params.aas} \
      --dels ${params.dels}
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
        uk_get_snps(uk_alignment.out)
        type_AAs_and_dels(uk_alignment.out, uk_metadata)
    emit:
        uk_get_variants.out
        uk_get_indels.out.uk_insertions
        uk_get_indels.out.uk_deletions
        fasta = uk_alignment.out
        metadata = type_AAs_and_dels.out
}


workflow {
    align_and_variant_call_cog_uk(params.uk_fasta,
                                  params.uk_metadata)
}