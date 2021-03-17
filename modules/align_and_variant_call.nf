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
    * Creates TSV of indels found in each genome
    * @input sam
    * @output insertions, deletions
    */

    publishDir "${publish_dev}/", pattern: "*/*.tsv", mode: 'copy'
    publishDir "${publish_dir}/", pattern: "*/*.tsv", mode: 'copy', enabled: { ${category} == 'cog'}

    input:
    path sam
    val category

    output:
    path "${category}/${category}.insertions.tsv", emit: insertions
    path "${category}/${category}.deletions.tsv", emit: deletions

    script:
    """
    mkdir -p ${category}
    gofasta sam indels \
      -s ${sam} \
      --threshold 2 \
      --insertions-out "${category}/${category}.insertions.tsv" \
      --deletions-out "${category}/${category}.deletions.tsv"
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
    mkdir -p ${category}
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

process get_nuc_variants {
    /**
    * Combines nucleotide variants into a metadata file which can be merged into the master
    * @input snps, dels
    * @output metadata
    */

    input:
    path snps
    path dels

    output:
    path "nuc_variants.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    sample_dict = {}
    with open("${dels}", 'r', newline = '') as csv_in:
        for line in csv_in:
            ref_start, length, samples = line.strip().split()
            samples = samples.split('|')
            var = "del_%s_%s" %(ref_start, length)
            for sample in samples:
                if sample in sample_dict:
                    sample_dict[sample].append(var)
                else:
                    sample_dict[sample] = [var]

    with open("${snps}", 'r', newline = '') as csv_in, \
        open("nuc_variants.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = ["sequence_name", "nucleotide_variants"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            row["sequence_name"] = row["query"]
            row["nucleotide_variants"] = row["SNPs"]
            if row["sequence_name"] in sample_dict:
                all_vars = [row["nucleotide_variants"]]
                all_vars.extend(sample_dict[row["sequence_name"]])
                row["nucleotide_variants"] = '|'.join(all_vars)
            for key in [k for k in row if k not in ["sequence_name", "nucleotide_variants"]]:
                del row[key]
            writer.writerow(row)
    """
}


process add_nucleotide_variants_to_metadata {
    /**
    * Adds nucleotide variants to metadata
    * @input metadata, nucleotide_variants
    * @output metadata
    */

    publishDir "${publish_dev}/cog_gisaid", pattern: "*.csv", mode: 'copy', saveAs: {"cog_gisaid_master.csv"}

    input:
    path metadata
    path nucleotide_variants

    output:
    path "${metadata.baseName}.with_nuc_variants.csv"

    script:
    """
    fastafunk add_columns \
          --in-metadata ${metadata} \
          --in-data ${nucleotide_variants} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns nucleotide_variants \
          --out-metadata "${metadata.baseName}.with_nuc_variants.csv"
    """
}

workflow align_and_variant_call {
    take:
        in_fasta
        in_metadata
        category
    main:
        minimap2_to_reference(in_fasta)
        get_variants(minimap2_to_reference.out, category)
        get_indels(minimap2_to_reference.out, category)
        alignment(minimap2_to_reference.out)
        mask_alignment(alignment.out)
        get_snps(mask_alignment.out, category)
        type_AAs_and_dels(mask_alignment.out, get_variants.out, category)
        get_nuc_variants(get_snps.out, get_indels.out.deletions)
        add_nucleotide_variants_to_metadata(in_metadata, get_nuc_variants.out)
    emit:
        variants = type_AAs_and_dels.out
        fasta = mask_alignment.out
        metadata = add_nucleotide_variants_to_metadata.out
}


mask_file = file(params.mask_file)
aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)


workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    category = params.category

    align_and_variant_call(uk_fasta, uk_metadata, category)
}