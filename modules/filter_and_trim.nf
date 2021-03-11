#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process filter_low_coverage_sequences {
    /**
    * Keeps only sequences with completeness greater than min_covg threshold
    * @input alignment, metadata
    * @output alignment_updated, metadata_updated
    * @params min_covg
    */

    input:
    path alignment
    path metadata

    output:
    path "${alignment.baseName}.low_covg_filtered.fasta", emit: fasta_updated
    path "${metadata.baseName}.low_covg_filtered.csv", emit: metadata_updated

    script:
    if (!params.min_covg)
        """
        mv "${alignment}" "${alignment.baseName}.low_covg_filtered.fasta"
        mv "${metadata}" "${metadata.baseName}.low_covg_filtered.csv"
        """
    else
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${alignment}", "fasta")

        with open("${metadata}", 'r', newline = '') as csv_in, \
             open("${metadata.baseName}.low_covg_filtered.csv", 'w', newline = '') as csv_out, \
             open("${alignment.baseName}.low_covg_filtered.fasta", 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row["why_excluded"]:
                    writer.writerow(row)
                    continue
                id = row["sequence_name"]
                if id in alignment:
                    seq = str(alignment[id].seq)
                    mapped_completeness = float(len(seq.replace("N", "")) / len(seq))
                    if mapped_completeness >= float(${params.min_covg} / 100):
                        writer.writerow(row)
                        fasta_out.write(">" + id + "\\n")
                        fasta_out.write(seq + "\\n")
                    else:
                        row["why_excluded"] = "low mapped_completeness"
                        writer.writerow(row)
        """
}


process publish_master_metadata {
    /**
    * Publishes master metadata csv for this category
    * @input metadata
    * @output metadata
    */

    publishDir "${publish_dev}", pattern: "*/*.csv", mode: 'copy'

    input:
    path metadata
    val category

    output:
    path "${category}/${category}_master.csv"

    script:
    """
    mkdir -p ${category}
    cp ${metadata} ${category}/${category}_master.csv
    """
}


process trim_alignment {
    /**
    * Trims start and end of alignment
    * @input alignment
    * @output alignment_updated
    * @params trim_start, trim_end
    */

    input:
    path alignment

    output:
    path "${alignment.baseName}.trimmed.fa"

    script:
    if (params.trim_start && params.trim_end)
        """
        #!/usr/bin/env python3
        from Bio import SeqIO

        strt = int(${params.trim_start})
        stp = int(${params.trim_end})

        with open("${alignment}", "r") as fasta_in, \
             open("${alignment.baseName}.trimmed.fa", "w") as fasta_out:

            for record in SeqIO.parse(fasta_in, "fasta"):
                seq = str(record.seq).upper()
                new_seq = ("N" * strt) + seq[strt:stp] + ("N" * (len(seq) - stp))
                fasta_out.write(">" + record.id + "\\n")
                fasta_out.write(new_seq + "\\n")
        """
    else
        """
        mv "${alignment.baseName}" "${alignment.baseName}.trimmed.fa"
        """
}


process distance_QC {
    /**
    * Outputs number of sequences per country
    * @input fasta, metadata
    * @output "QC_distances.tsv"
    */
    publishDir "${publish_dev}", pattern: "*/*.tsv", mode: 'copy'


    input:
    path fasta
    path metadata
    val category

    output:
    path "${category}/${category}_QC_distances.tsv"

    script:
    """
    datafunk distance_to_root \
        --input-fasta ${fasta} \
        --input-metadata ${metadata}

    mkdir -p ${category}
    mv distances.tsv "${category}/${category}_QC_distances.tsv"
    """
}


process filter_on_distance_to_WH04 {
    /**
    * Restricts to samples within distance x of WH04
    * @input fasta, metadata, distances
    * @output
    */

    input:
    path fasta
    path metadata
    path distances

    output:
    path "${fasta.baseName}.distance_filtered.fa", emit: fasta
    path "${metadata.baseName}.distance_filtered.csv", emit: metadata

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    reject = []
    with open("${distances}", 'r', newline = '') as distances_in:
        reader = csv.DictReader(distances_in, delimiter="\t", quotechar='\"', dialect = "unix")
        for row in reader:
            sequence_name = row['sequence_name']
            distance = float(row['distance_stdevs'])
            if distance >= 4.0:
                reject.append(sequence_name)

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${metadata.baseName}.distance_filtered.csv", 'w', newline = '') as csv_out, \
        open("${fasta.baseName}.distance_filtered.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            if row["why_excluded"]:
                writer.writerow(row)
                continue
            id = row["sequence_name"]
            if id in reject:
                row["why_excluded"] = "distance to WH04 more than 4.0 epi-week std devs"
                writer.writerow(row)
                continue
            if id in alignment:
                writer.writerow(row)
                seq = str(alignment[id].seq)
                fasta_out.write(">" + id + "\\n")
                fasta_out.write(seq + "\\n")
    """
}


workflow filter_and_trim_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
    main:
        filter_low_coverage_sequences(gisaid_fasta, gisaid_metadata)
        trim_alignment(filter_low_coverage_sequences.out.fasta_updated)
        distance_QC(trim_alignment.out, filter_low_coverage_sequences.out.metadata_updated, "gisaid")
        filter_on_distance_to_WH04(gisaid_fasta, gisaid_metadata, distance_QC.out)
        publish_master_metadata(filter_on_distance_to_WH04.out.metadata, "gisaid")
    emit:
        fasta = filter_on_distance_to_WH04.out.fasta
        metadata = filter_on_distance_to_WH04.out.metadata
}


workflow filter_and_trim_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        filter_low_coverage_sequences(uk_fasta, uk_metadata)
        trim_alignment(filter_low_coverage_sequences.out.fasta_updated)
        publish_master_metadata(filter_low_coverage_sequences.out.metadata_updated, "cog")
    emit:
        fasta = trim_alignment.out
        metadata = filter_low_coverage_sequences.out.metadata_updated
}

workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)

    filter_and_trim_cog_uk(uk_fasta,
                           uk_metadata)
}
