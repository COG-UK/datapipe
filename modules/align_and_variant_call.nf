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
    minimap2 -t ${task.cpus} -a --secondary=no -x asm20 --score-N=0 ${reference_fasta} ${fasta} > alignment.sam
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
    label 'retry_increasing_mem'


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

process get_updown {
    /**
    * Call SNPs in each genome
    * @input alignment
    * @output updown list
    * @params reference_fasta
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path alignment
    val category

    output:
    path "${category}/${category}.updown.csv"

    script:
    """
    mkdir -p ${category}
    gofasta updown list -r ${WH04_fasta} -q ${alignment} -o ${category}/${category}.updown.csv
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
    sed "s/query/sequence_name/g" "mutations.tmp.csv" > mutations.tmp2.csv
    sed "s/variants/mutations/g" "mutations.tmp2.csv" > ${category}/${category}_mutations.csv

    if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${category}/${category}_mutations.csv" | wc -l) ]]
        then
            echo \$(cat "${metadata}" | wc -l)
            echo \$(cat "${category}/${category}_mutations.csv" | wc -l)
            exit 1
        fi
    """
}

process get_nuc_mutations {
    /**
    * Combines nucleotide mutations into a metadata file which can be merged into the master
    * @input snps, dels, ins
    * @output metadata
    */

    input:
    path snps
    path dels
    path ins

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

    with open("${ins}", 'r', newline = '') as csv_in:
        for line in csv_in:
            ref_start, insertion, samples= line.strip().split()
            samples = samples.split('|')
            var = "ins_%s_%s" %(ref_start, insertion)
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


process restrict_metadata {
    /**
    * restricts only to sequences not excluded
    * @input metadata
    * @output metadata
    */

    input:
    path metadata

    output:
    path "${metadata.baseName}.restricts.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${metadata.baseName}.restricts.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            if row["why_excluded"] not in [None, "", "None"]:
                writer.writerow(row)
    """
}


process add_nucleotide_mutations_to_metadata {
    /**
    * Adds nucleotide mutations to metadata
    * @input metadata, nucleotide_mutations
    * @output metadata
    */

    memory { 1.GB * task.attempt + metadata.size() * 2.B }

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

    if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${metadata.baseName}.with_nuc_mutations.csv" | wc -l) ]]
        then
            echo \$(cat "${metadata}" | wc -l)
            echo \$(cat "${metadata.baseName}.with_nuc_mutations.csv" | wc -l)
            exit 1
        fi
    """
}

process add_ambiguities_to_metadata {
    /**
    * Adds nucleotide mutations to metadata
    * @input metadata, nucleotide_mutations
    * @output metadata
    */

    memory { 1.GB * task.attempt + metadata.size() * 2.B }

    input:
    path metadata
    path updown

    output:
    path "${metadata.baseName}.with_ambiguities.csv"

    script:
    """
    fastafunk add_columns \
          --in-metadata ${metadata} \
          --in-data ${updown} \
          --index-column sequence_name \
          --join-on querye \
          --new-columns ambiguities \
          --out-metadata "${metadata.baseName}.with_ambiguities.csv"

    if [[ \$(cat "${metadata}" | wc -l) != \$(cat "${metadata.baseName}.with_ambiguities.csv" | wc -l) ]]
        then
            echo \$(cat "${metadata}" | wc -l)
            echo \$(cat "${metadata.baseName}.with_ambiguities.csv" | wc -l)
            exit 1
        fi
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
      -n ${params.constellations}

    if [[ \$(grep ">" "${alignment}" | wc -l) != \$(tail -n+2 "${alignment.baseName}.haplotyped.csv" | wc -l) ]]
            then
                echo \$(grep ">" "${alignment}" | wc -l)
                echo \$(tail -n+2 "${alignment.baseName}.haplotyped.csv" | wc -l)
                exit 1
            fi
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
      -n ${params.constellations}

    if [[ \$(grep ">" "${alignment}" | wc -l) != \$(tail -n+2 "${alignment.baseName}.classified.csv" | wc -l) ]]
                then
                    echo \$(grep ">" "${alignment}" | wc -l)
                    echo \$(tail -n+2 "${alignment.baseName}.classified.csv" | wc -l)
                    exit 1
                fi
    """
}

process add_constellations_to_metadata {
    /**
    * Adds constellations to metadata
    * @input metadata, haplotyped, classified
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'
 
    memory { task.attempt * (classified.size() + haplotyped.size()) * 9.B }

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

    if [[ \$(cat "${haplotyped}" | wc -l) != \$(cat "${category}/${category}_constellations.csv" | wc -l) ]]
            then
                echo \$(cat "${haplotyped}" | wc -l)
                echo \$(cat "${category}/${category}_constellations.csv" | wc -l)
                exit 1
            fi
    """
}


process announce_summary {
    /**
    * Summarizes alignment into JSON
    * @input fastas
    */

    input:
    path fasta
    path alignment

    output:
    path "announce.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Finished alignment and variant calling ${params.date}*\\n" >> announce.json
                echo "> Number of sequences in FASTA : \$(cat ${fasta} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences in ALIGNMENT : \$(cat ${alignment} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Finished alignment and variant calling ${params.date}*\\n" >> announce.json
                echo "> Number of sequences in FASTA : \$(cat ${fasta} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences in ALIGNMENT : \$(cat ${alignment} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json
            """
}

workflow align_and_variant_call {
    take:
        in_fasta
        in_metadata
        category
    main:
        in_fasta.splitFasta( by: params.chunk_size, file: true ).set{ fasta_chunks }
        minimap2_to_reference(fasta_chunks)
        alignment(minimap2_to_reference.out)
        alignment.out.collectFile(newLine: false).set{ alignment_result }
        minimap2_to_reference.out.collectFile(newLine: false, keepHeader: true, skip: 2).set{ mapped_result }


        get_mutations(mapped_result, category)
        get_indels(mapped_result, category)

        get_snps(alignment_result, category)
        get_updown(alignment_result, category)
        type_AAs_and_dels(alignment_result, get_mutations.out, category)
        get_nuc_mutations(get_snps.out, get_indels.out.deletions, get_indels.out.insertions)
        add_nucleotide_mutations_to_metadata(in_metadata, get_nuc_mutations.out)
        add_ambiguities_to_metadata(type_AAs_and_dels.out, get_updown.out)

        haplotype_constellations(alignment.out)
        haplotype_constellations.out.collectFile(newLine: false, keepHeader: true, skip: 1).set{ haplotype_result }
        classify_constellations(alignment.out)
        classify_constellations.out.collectFile(newLine: false, keepHeader: true, skip: 1).set{ classify_result }

        add_constellations_to_metadata(haplotype_result, classify_result, category)
        announce_summary(in_fasta, alignment_result)
    emit:
        mutations = add_ambiguities_to_metadata.out
        constellations = add_constellations_to_metadata.out
        fasta = alignment_result
        metadata = add_nucleotide_mutations_to_metadata.out
        updown = get_updown.out
}


aas = file(params.aas)
dels = file(params.dels)
reference_fasta = file(params.reference_fasta)
reference_genbank = file(params.reference_genbank)
WH04_fasta = file(params.WH04_fasta)

workflow {
    uk_fasta = Channel.fromPath(params.uk_fasta)
    uk_metadata = Channel.fromPath(params.uk_metadata)
    category = params.category

    align_and_variant_call(uk_fasta, uk_metadata, category)
}
