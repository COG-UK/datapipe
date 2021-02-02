#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dir = file(params.publish_dir)


process uk_filter_low_coverage_sequences {
    /**
    * Keeps only sequences with completeness greater than min_covg threshold
    * @input uk_alignment, uk_metadata
    * @output uk_alignment_updated, uk_metadata_updated
    * @params min_covg
    */

    input:
    file uk_alignment
    file uk_metadata

    output:
    path "${uk_alignment.baseName}.low_covg_filtered.fasta", emit: uk_fasta_updated
    path "${uk_metadata.baseName}.low_covg_filtered.csv", emit: uk_metadata_updated

    script:
    if (!params.min_covg)
        """
        mv "${uk_alignment}" "${uk_alignment.baseName}.low_covg_filtered.fasta"
        mv "${uk_metadata}" "${uk_metadata.baseName}.low_covg_filtered.csv"
        """
    else
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${uk_alignment}", "fasta")

        with open("${uk_metadata}", 'r', newline = '') as csv_in, \
             open("${uk_metadata.baseName}.low_covg_filtered.csv", 'w', newline = '') as csv_out, \
             open("${uk_alignment.baseName}.low_covg_filtered.fasta", 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                id = row["fasta_header"]
                if id in alignment:
                    seq = str(alignment[id].seq)
                    mapped_completeness = float(len(seq.replace("N", "")) / len(seq))
                    if mapped_completeness >= float(${params.min_covg} / 100):
                        writer.writerow(row)
                        fasta_out.write(">" + id + "\\n")
                        fasta_out.write(seq + "\\n")
                    else:
                        continue
        """
}


process uk_trim_alignment {
    /**
    * Trims start and end of alignment
    * @input uk_alignment
    * @output uk_alignment_updated
    * @params trim_start, trim_end
    */

    input:
    file uk_alignment

    output:
    path "${uk_alignment.baseName}.trimmed.fa"

    script:
    if (params.trim_start && params.trim_end)
        """
        #!/usr/bin/env python3
        from Bio import SeqIO

        strt = int(${params.trim_start})
        stp = int(${params.trim_end})

        with open("${uk_alignment}", "r") as fasta_in, \
             open("${uk_alignment.baseName}.trimmed.fa", "w") as fasta_out:

            for record in SeqIO.parse(fasta_in, "fasta"):
                seq = str(record.seq).upper()
                new_seq = ("N" * strt) + seq[strt:stp] + ("N" * (len(seq) - stp))
                fasta_out.write(">" + record.id + "\\n")
                fasta_out.write(new_seq + "\\n")
        """
    else
        """
        mv "${uk_alignment.baseName}" "${uk_alignment.baseName}.trimmed.fa"
        """
}

process publish_filtered_aligned_cog_data {
    /**
    * Publish filtered alignment and relevant metadata
    * @input uk_alignment, uk_metadata
    * @params date
    */

    publishDir "${publish_dir}/alignments/", pattern: "*.fa", mode: 'copy', saveAs: {"cog_${params.date}_alignment.fasta"}
    publishDir "${publish_dir}/alignments/", pattern: "*.csv", mode: 'copy', saveAs: {"cog_${params.date}_metadata.fasta"}

    input:
    file uk_alignment
    file uk_metadata

    output:
    path "${uk_alignment.baseName}.matched.fa"
    path "${uk_metadata.baseName}.matched.csv"

    script:
    """
    fastafunk fetch \
          --in-fasta ${uk_alignment} \
          --in-metadata ${uk_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name secondary_identifier sample_date epi_week \
                                    country adm1 adm2 outer_postcode \
                                    is_surveillance is_community is_hcw \
                                    is_travel_history travel_history lineage \
                                    lineage_support \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private \
          --out-fasta "${uk_alignment.baseName}.matched.fa" \
          --out-metadata "${uk_metadata.baseName}.matched.csv" \
          --restrict
    """
}

workflow filter_and_trim_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        uk_filter_low_coverage_sequences(uk_fasta, uk_metadata)
        uk_trim_alignment(uk_filter_low_coverage_sequences.out.uk_fasta_updated)
        publish_filtered_aligned_cog_data(uk_trim_alignment.out,uk_filter_low_coverage_sequences.out.uk_metadata_updated)
    emit:
        fasta = uk_trim_alignment.out
        metadata = uk_filter_low_coverage_sequences.out.uk_metadata_updated
}

workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)

    filter_and_trim_cog_uk(uk_fasta,
                           uk_metadata)
}
