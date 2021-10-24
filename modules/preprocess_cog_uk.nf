#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process uk_strip_header_digits_and_unalign {
    /**
    * Strips extra header info from FASTA, removed '-' from sequence
    * @input uk_fasta
    * @output uk_fasta_updated
    */

    input:
    path uk_fasta

    output:
    path "${uk_fasta.baseName}.header_stripped.fasta"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    fasta_in = SeqIO.parse("${uk_fasta}", "fasta")
    with open("${uk_fasta.baseName}.header_stripped.fasta", 'w') as f:
        for record in fasta_in:
            ID = record.description.split("|")[0]
            f.write(">" + ID + "\\n")
            seq = str(record.seq).replace('-','')
            seq = seq.replace('?','N')
            f.write(seq + "\\n")
    """
}


process uk_add_columns_to_metadata {
    /**
    * Takes the MAJORA TSV of metadata and adds/updates columns for sample_date, pillar_2,
    * sequence_name, covv_accession_id, edin_epi_week, edin_epi_day and adm0
    * @input uk_metadata
    * @output uk_metadata_updated
    * @params uk_accessions, uk_updated_dates
    */

    input:
    path uk_metadata
    path uk_accessions
    path uk_updated_dates

    output:
    path "${uk_metadata.baseName}.updated.csv"

    script:
    """
    $project_dir/../bin/add_to_uk_metadata.py \
        --in-metadata ${uk_metadata} \
        --out-metadata ${uk_metadata.baseName}.updated.csv \
        --accession-file ${uk_accessions} \
        --updated-date-file ${uk_updated_dates}
    """
}

process uk_filter_omitted_sequences {
    /**
    * Takes a FASTA and METADATA and excludes samples specified in an exclusion file
    * sequence_name, covv_accession_id, edin_epi_week, edin_epi_day and adm0
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    * @params uk_omissions
    */
    input:
    path uk_fasta
    path uk_metadata
    path uk_omissions

    output:
    path "${uk_fasta.baseName}.omit_filtered.fa", emit: fasta
    path "${uk_metadata.baseName}.omit_filtered.csv", emit: metadata

    script:
    if ( params.uk_omissions )
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${uk_fasta}", "fasta")

        omissions = set()
        with open("${uk_omissions}", "r") as f:
            for line in f:
                omissions.add(line.rstrip())

        with open("${uk_metadata}", 'r', newline = '') as csv_in, \
             open("${uk_metadata.baseName}.omit_filtered.csv", 'w', newline = '') as csv_out, \
             open("${uk_fasta.baseName}.omit_filtered.fa", "w") as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row["central_sample_id"] in omissions:
                    row["why_excluded"] = "central_sample_id in omissions_file"
                    writer.writerow(row)
                    continue

                record = alignment[row["fasta_header"]]
                writer.writerow(row)
                fasta_out.write(">" + record.id + "\\n")
                fasta_out.write(str(record.seq) + "\\n")
        """
    else
        """
        mv "${uk_fasta}" "${uk_fasta.baseName}.omit_filtered.fa"
        mv "${uk_metadata}" "${uk_metadata.baseName}.omit_filtered.csv"
        """
}

process uk_filter_on_sample_date {
    /**
    * If a time window (in days) is provided, excludes samples from FASTA and
    * METADATA files which do not fall within X days of date
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_update, uk_metadata_updated
    * @params time_window, date
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.date_filtered.fa", emit: fasta
    path "${uk_metadata.baseName}.date_filtered.csv", emit: metadata

    script:
    if ( params.time_window && params.date)
        """
        #!/usr/bin/env python3
        import datetime
        from Bio import SeqIO
        import csv

        indexed_fasta = SeqIO.index("${uk_fasta}", "fasta")

        window = datetime.timedelta(int("${params.time_window}"))
        todays_date = datetime.datetime.strptime("${params.date}", '%Y-%m-%d').date()

        with open"${uk_metadata}", 'r', newline = '') as csv_in, \
            open("${uk_metadata.baseName}.date_filtered.csv", 'w', newline = '') as csv_out, \
            open("${uk_fasta.baseName}.date_filtered.fa", "w") as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                try:
                    date = datetime.datetime.strptime(row["sample_date"], '%Y-%m-%d').date()
                except:
                    row["why_excluded"] = "no sample_date"
                    writer.writerow(row)

                if (todays_date - window) > date:
                    row["why_excluded"] = "sample_date older than %s days" %window
                    writer.writerow(row)
                    continue

                writer.writerow(row)

                seq_rec = indexed_fasta[row["fasta_header"]]
                fasta_out.write(">" + seq_rec.id + "\\n")
                fasta_out.write(str(seq_rec.seq) + "\\n")
        """
    else
        """
        mv "${uk_fasta}" "${uk_fasta.baseName}.date_filtered.fa"
        mv "${uk_metadata}" "${uk_metadata.baseName}.date_filtered.csv"
        """
}


process add_previous_uk_lineage_to_metadata {
    /**
    * Adds uk_lineage where previously assigned
    * @input metadata
    * @output metadata
    */
    label 'retry_increasing_mem'

    input:
    path metadata

    output:
    path "${metadata.baseName}.with_uk_lineage.csv"

    script:
    if ( !params.previous_metadata )
        """
        mv ${metadata} "${metadata.baseName}.with_uk_lineage.csv"
        """
    else
        """
        fastafunk add_columns \
              --in-metadata ${metadata} \
              --in-data ${params.previous_metadata} \
              --index-column sequence_name \
              --join-on sequence_name \
              --new-columns uk_lineage \
              --out-metadata "${metadata.baseName}.with_uk_lineage.csv"
        """
}


process announce_summary {
    /**
    * Summarizes preprocess into JSON
    * @input fastas
    */

    input:
    path original
    path strip_header
    path filter_omitted_sequences
    path filter_on_sample_date

    output:
    path "announce.json"

    script:
        if (params.webhook)
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Preprocessing COG input ${params.date}*\\n" >> announce.json
                echo "> Number of sequences in COG input files : \$(cat ${original} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after header stripped : \$(cat ${strip_header} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering omitted: \$(cat ${filter_omitted_sequences} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering by sample date with time window ${params.time_window}: \$(cat ${filter_on_sample_date} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json

            echo 'webhook ${params.webhook}'

            curl -X POST -H "Content-type: application/json" -d @announce.json ${params.webhook}
            """
        else
            """
            echo '{"text":"' > announce.json
                echo "*${params.whoami}: Preprocessing COG input ${params.date}*\\n" >> announce.json
                echo "> Number of sequences in COG input files : \$(cat ${original} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after header stripped : \$(cat ${strip_header} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering omitted: \$(cat ${filter_omitted_sequences} | grep '>' | wc -l)\\n" >> announce.json
                echo "> Number of sequences after filtering by sample date with time window ${params.time_window}: \$(cat ${filter_on_sample_date} | grep '>' | wc -l)\\n" >> announce.json
                echo '"}' >> announce.json
            """
}

uk_updated_dates = file(params.uk_updated_dates)
uk_omissions = file(params.uk_omissions)

workflow preprocess_cog_uk {
    take:
        uk_fasta
        uk_metadata
        uk_accessions
    main:
        uk_strip_header_digits_and_unalign(uk_fasta)
        uk_add_columns_to_metadata(uk_metadata, uk_accessions, uk_updated_dates)
        uk_filter_omitted_sequences(uk_strip_header_digits_and_unalign.out, uk_add_columns_to_metadata.out, uk_omissions)
        uk_filter_on_sample_date(uk_filter_omitted_sequences.out.fasta, uk_filter_omitted_sequences.out.metadata)
        add_previous_uk_lineage_to_metadata(uk_filter_omitted_sequences.out.metadata)
        announce_summary(uk_fasta, uk_strip_header_digits_and_unalign.out, uk_filter_omitted_sequences.out.fasta, uk_filter_on_sample_date.out.fasta)
    emit:
        fasta = uk_filter_on_sample_date.out.fasta
        metadata = add_previous_uk_lineage_to_metadata.out
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    uk_accessions = file(params.uk_accessions)

    preprocess_cog_uk(uk_fasta,
                      uk_metadata,
                      uk_accessions)
}
