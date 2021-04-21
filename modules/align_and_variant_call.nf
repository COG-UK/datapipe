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

process get_mutations {
    /**
    * Creates CSV of mutations found in each genome
    * @input sam
    * @output mutations
    * @parms reference_fasta, reference_genbank
    */

    cpus 4

    input:
    path sam
    val category

    output:
    path "${category}.mutations.csv"

    script:
    """
    gofasta sam variants -t ${task.cpus} \
      --samfile ${sam} \
      --reference ${reference_fasta} \
      --genbank ${reference_genbank} \
      --outfile ${category}.mutations.csv
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
    path "${category}/${category}_mutations.csv"

    script:
    """
    mkdir -p ${category}
    $project_dir/../bin/type_aas_and_dels.py \
      --in-fasta ${alignment} \
      --in-metadata ${metadata} \
      --out-metadata "mutations.tmp.csv" \
      --reference-fasta ${reference_fasta} \
      --aas ${aas} \
      --dels ${dels} \
      --index-column query
    sed "s/query/sequence_name/g" "mutations.tmp.csv" > ${category}/${category}_mutations.csv

    """
}

process get_nuc_mutations {
    /**
    * Combines nucleotide mutations into a metadata file which can be merged into the master
    * @input snps, dels
    * @output metadata
    */

    input:
    path snps
    path dels

    output:
    path "nuc_mutations.csv"

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
        open("nuc_mutations.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = ["sequence_name", "nucleotide_mutations"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            row["sequence_name"] = row["query"]
            row["nucleotide_mutations"] = row["SNPs"]
            if row["sequence_name"] in sample_dict:
                all_vars = [row["nucleotide_mutations"]]
                all_vars.extend(sample_dict[row["sequence_name"]])
                row["nucleotide_mutations"] = '|'.join(all_vars)
            for key in [k for k in row if k not in ["sequence_name", "nucleotide_mutations"]]:
                del row[key]
            writer.writerow(row)
    """
}


process add_nucleotide_mutations_to_metadata {
    /**
    * Adds nucleotide mutations to metadata
    * @input metadata, nucleotide_mutations
    * @output metadata
    */

    input:
    path metadata
    path nucleotide_mutations

    output:
    path "${metadata.baseName}.with_nuc_mutations.csv"

    script:
    """
    fastafunk add_columns \
          --in-metadata ${metadata} \
          --in-data ${nucleotide_mutations} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns nucleotide_mutations \
          --out-metadata "${metadata.baseName}.with_nuc_mutations.csv"
    """
}


process haplotype_constellations {
    /**
    * Adds a column to metadata table for each constellation, and a summary column for all found
    * @input alignment
    * @output haplotype_csv
    * @params constellations
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.haplotyped.csv"

    script:
    """
    scorpio haplotype \
      --input ${alignment} \
      --output "${alignment.baseName}.haplotyped.csv" \
      --output-counts \
      --constellations ${constellations}/*.json
    """
}

process classify_constellations {
    /**
    * Adds a column to metadata table for each constellation, and a summary column for all found
    * @input alignment
    * @output classify_csv
    * @params constellations
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.classified.csv"

    script:
    """
    scorpio classify \
      --input ${alignment} \
      --output "${alignment.baseName}.classified.csv" \
      --constellations ${constellations}/*.json
    """
}

process add_constellations_to_metadata {
    /**
    * Adds constellations to metadata
    * @input metadata, haplotyped, classified
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path haplotyped
    path classified
    val category

    output:
    path "${category}/${category}_constellations.csv"

    script:
    """
    mkdir -p ${category}
    fastafunk add_columns \
          --in-metadata ${classified} \
          --in-data ${haplotyped} \
          --index-column query \
          --join-on query \
          --out-metadata "constellations.tmp.csv"
    sed "s/query/sequence_name/g" "constellations.tmp.csv" > "${category}/${category}_constellations.csv"
    """
}

workflow align_and_variant_call {
    take:
        in_fasta
        in_metadata
        category
    main:
        minimap2_to_reference(in_fasta)
        get_mutations(minimap2_to_reference.out, category)
        get_indels(minimap2_to_reference.out, category)
        alignment(minimap2_to_reference.out)
        get_snps(alignment.out, category)
        type_AAs_and_dels(alignment.out, get_mutations.out, category)
        get_nuc_mutations(get_snps.out, get_indels.out.deletions)
        add_nucleotide_mutations_to_metadata(in_metadata, get_nuc_mutations.out)
        haplotype_constellations(alignment.out)
        classify_constellations(alignment.out)
        add_constellations_to_metadata(haplotype_constellations.out, classify_constellations.out, category)
    emit:
        mutations = type_AAs_and_dels.out
        constellations = add_constellations_to_metadata.out
        fasta = alignment.out
        metadata = add_nucleotide_mutations_to_metadata.out
}


aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)
constellations = file(params.constellations)

workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    category = params.category

    align_and_variant_call(uk_fasta, uk_metadata, category)
}