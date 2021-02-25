#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)
publish_dev = file(params.publish_dev)


process minimap2_to_reference {
    /**
    * Minimaps samples to reference
    * @input fasta
    * @output sam
    * @params reference_fasta
    */

    cpus 4

    input:
    path fasta

    output:
    path "alignment.sam"

    script:
    """
    minimap2 -t ${task.cpus} -a -x asm5 ${reference_fasta} ${fasta} > alignment.sam
    """
}

process get_variants {
    /**
    * Creates CSV of variants found in each genome
    * @input sam
    * @output variants
    * @parms reference_fasta, reference_genbank
    */

    cpus 4

    input:
    path sam
    val category

    output:
    path "${category}.variants.csv"

    script:
    """
    gofasta sam variants -t ${task.cpus} \
      --samfile ${sam} \
      --reference ${reference_fasta} \
      --genbank ${reference_genbank} \
      --outfile ${category}.variants.csv
    """
}

process get_indels {
    /**
    * Creates CSV of indels found in each genome
    * @input sam
    * @output insertions, deletions
    */

    publishDir "${publish_dev}/metadata/", pattern: "*.csv", mode: 'copy'

    input:
    path sam
    val category

    output:
    path "${category}.insertions.csv", emit: insertions
    path "${category}.deletions.csv", emit: deletions

    script:
    """
    gofasta sam indels \
      -s ${sam} \
      --threshold 2 \
      --insertions-out "${category}.insertions.csv" \
      --deletions-out "${category}.deletions.csv"
    """
}

process alignment {
    /**
    * Get reference-based alignment
    * @input sam
    * @output alignment
    * @params reference_fasta
    */

    cpus 4

    input:
    path sam

    output:
    path "alignment.fasta"

    script:
    """
    gofasta sam toMultiAlign -t ${task.cpus} \
      --samfile ${sam} \
      --reference ${reference_fasta} \
      --pad \
      -o alignment.fasta
    """
}


process mask_alignment {
    /**
    * Applies a mask to aligned FASTA
    * @input alignment
    * @output alignment_updated
    * @params mask_file
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.masked.fa"

    script:
    """
    $project_dir/../bin/add_mask.py \
      --in-alignment ${alignment} \
      --out-alignment "${alignment.baseName}.masked.fa" \
      --mask ${mask_file} \
    """
}


process get_snps {
    /**
    * Call SNPs in each genome
    * @input alignment
    * @output snps
    * @params reference_fasta
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path alignment
    val category

    output:
    path "${category}/${category}.snps.csv"

    script:
    """
    mkdir ${category}
    gofasta snps -r ${reference_fasta} -q ${alignment} -o ${category}/${category}.snps.csv
    """
}

process type_AAs_and_dels {
    /**
    * Adds a column to metadata table for specific dels and aas looked for
    * @input alignment, metadata
    * @output metadata_updated
    * @params reference_fasta, del_file, aa_file
    */

    publishDir "${publish_dev}/", pattern: "*/*.csv", mode: 'copy'

    input:
    path alignment
    path metadata
    val category

    output:
    path "${category}/${category}_variants.csv"

    script:
    """
    mkdir -p ${category}
    $project_dir/../bin/type_aas_and_dels.py \
      --in-fasta ${alignment} \
      --in-metadata ${metadata} \
      --out-metadata ${category}/${category}_variants.csv \
      --reference-fasta ${reference_fasta} \
      --aas ${aas} \
      --dels ${dels} \
      --index-column query
    """
}

workflow align_and_variant_call {
    take:
        in_fasta
        category
    main:
        minimap2_to_reference(in_fasta)
        get_variants(minimap2_to_reference.out, category)
        get_indels(minimap2_to_reference.out, category)
        alignment(minimap2_to_reference.out)
        mask_alignment(alignment.out)
        get_snps(mask_alignment.out, category)
        type_AAs_and_dels(mask_alignment.out, get_variants.out, category)
    emit:
        variants = type_AAs_and_dels.out
        fasta = mask_alignment.out
}


mask_file = file(params.mask_file)
aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)


workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    category = params.category

    align_and_variant_call_cog_uk(uk_fasta, category)
}