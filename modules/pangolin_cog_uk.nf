#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process extract_sequences_for_pangolin {
    /**
    * If update_all_lineage_assignments flag set, or no previous provided, outputs the input files.
    * Otherwise, extracts lineageless sequences from FASTA to run pangolin on, and updates
    * metadata with previous lineages
    * @input uk_fasta, uk_metadata
    * @output pangolin_fasta, metadata_with_previous
    * @params uk_previous_metadata, update_all_lineage_assignments
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.for_pangolin.fa", emit: pangolin_fasta
    path "${uk_metadata.baseName}.with_previous.csv", emit: metadata_with_previous

    script:
    if (params.update_all_lineage_assignments || !params.uk_previous_metadata )
        """
        mv "${uk_fasta}" "${uk_fasta.baseName}.for_pangolin.fa"
        mv "${uk_metadata}" "${uk_metadata.baseName}.with_previous.csv"
        """
    else
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${uk_fasta}", "fasta")

        lineage_dict = {}
        with open("${params.uk_previous_metadata}", 'r', newline = '') as lineages_in:
            reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")
            for row in reader:
                #if row["taxon"] in lineage_dict:
                if row["fasta_header"] in lineage_dict:
                    print("%s occurs more than once in lineages input file" % row["taxon"])
                    continue
                #lineage_dict[row["taxon"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["pangoLEARN_version"], "probability": row["probability"]}
                lineage_dict[row["fasta_header"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["lineages_version"], "probability": row["lineage_support"]}

        with open("${uk_metadata}", 'r', newline = '') as csv_in, \
            open("${uk_metadata.baseName}.with_previous.csv", 'w', newline = '') as csv_out, \
            open("${uk_fasta.baseName}.for_pangolin.fa", 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["lineage", "pangoLEARN_version", "probability"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                lineage = ""
                pangoLEARN_version = ""
                probability = ""

                fasta_header = row["fasta_header"]

                if fasta_header in lineage_dict:
                    lineage = lineage_dict[fasta_header]["lineage"]
                    pangoLEARN_version = lineage_dict[fasta_header]["pangoLEARN_version"]
                    probability = lineage_dict[fasta_header]["probability"]
                else:
                    seqrec = alignment[fasta_header]
                    fasta_out.write(">" + seqrec.id + "\\n")
                    fasta_out.write(str(seqrec.seq) + "\\n")

                row["lineage"] = lineage
                row["pangoLEARN_version"] = pangoLEARN_version
                row["probability"] = probability

                writer.writerow(row)
        """
}

process uk_pangolin {
    /**
    * Runs PANGOLIN on input fasta
    * @input uk_fasta
    * @output pangolin_fasta
    */

    input:
    path uk_fasta

    output:
    path "pangolin/lineage_report.csv"

    script:
    """
    pangolin "${uk_fasta}" \
        --outdir pangolin \
        --tempdir pangolin_tmp
    """
}

process uk_add_new_pangolin_lineages_to_metadata {
    /**
    * Updates metadata with new PANGOLIN lineage assignments
    * @input uk_metadata, pangolin_csv
    * @output uk_metadata_updated
    */

    input:
    path uk_metadata
    path pangolin_csv

    output:
    path "${uk_metadata.baseName}.with_pangolin.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv

    lineage_dict = {}
    with open("${pangolin_csv}", 'r', newline = '') as lineages_in:
        reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")
        for row in reader:
            if row["taxon"] in lineage_dict:
                print("%s occurs more than once in lineages input file" % row["taxon"])
                continue
            lineage_dict[row["taxon"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["pangoLEARN_version"], "probability": row["probability"]}


    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.with_pangolin.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        column_names = reader.fieldnames + [col for col in ["lineage", "pangoLEARN_version", "probability"] if col not in reader.fieldnames]
        writer = csv.DictWriter(csv_out, fieldnames = column_names, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["fasta_header"]
            if fasta_header in lineage_dict:
                row["lineage"] = lineage_dict[fasta_header]["lineage"]
                row["pangoLEARN_version"] = lineage_dict[fasta_header]["pangoLEARN_version"]
                row["probability"] = lineage_dict[fasta_header]["probability"]
                writer.writerow(row)
    """
}

workflow pangolin_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        extract_sequences_for_pangolin(uk_fasta, uk_metadata)
        uk_pangolin(extract_sequences_for_pangolin.out.pangolin_fasta)
        uk_add_new_pangolin_lineages_to_metadata(extract_sequences_for_pangolin.out.metadata_with_previous, uk_pangolin.out)
    emit:
        metadata = uk_add_new_pangolin_lineages_to_metadata.out
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)

    pangolin_cog_uk(uk_fasta,
                    uk_metadata)
}
