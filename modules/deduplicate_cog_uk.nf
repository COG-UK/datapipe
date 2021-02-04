#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process uk_annotate_with_unmapped_genome_completeness {
    /**
    * Adds a column to metadata with proportion of genome which is complete
    * @input uk_fasta, uk_metadata
    * @output uk_metadata_updated
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_metadata.baseName}.annotated.csv"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${uk_fasta}", "fasta")

    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.annotated.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["unmapped_genome_completeness"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            id = row["fasta_header"]
            if id in alignment:
                seq = str(alignment[id].seq)
                completeness = float(len(seq.replace("N", "")) / len(seq))
                row["unmapped_genome_completeness"] = completeness
                writer.writerow(row)
            else:
                continue
    """
}

process uk_remove_duplicates_COGID_by_proportionN {
    /**
    * Where duplicate COGID, keeps the most complete
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.deduplicated_by_cogid.fa", emit: uk_fasta_updated
    path "${uk_metadata.baseName}.deduplicated_by_cogid.csv", emit: uk_metadata_updated

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${uk_fasta}", "fasta")

    dup_dict = {}
    tokeep = set()

    with open("${uk_metadata}", 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

        for row in reader:
            fasta_header = row["fasta_header"]
            id = row["central_sample_id"]
            completeness = float(row["unmapped_genome_completeness"])

            if id in dup_dict:
                if completeness > dup_dict[id]["completeness"]:
                    dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}
                else:
                    continue
            else:
                dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}

    for k,v in dup_dict.items():
        tokeep.add(v["fasta_header"])

    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.deduplicated_by_cogid.csv", 'w', newline = '') as csv_out, \
         open("${uk_fasta.baseName}.deduplicated_by_cogid.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["fasta_header"]

            if fasta_header in tokeep:
                writer.writerow(row)
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\\n")
                fasta_out.write(str(seqrec.seq) + "\\n")
            else:
                continue
    """
}

process uk_remove_duplicates_biosamplesourceid_by_date {
    /**
    * Where duplicate biosample_source_id, keeps the earliest
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.deduplicated_by_biosamplesourceid.fa", emit: uk_fasta_updated
    path "${uk_metadata.baseName}.deduplicated_by_biosamplesourceid.csv", emit: uk_metadata_updated

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${uk_fasta}", "fasta")

    dup_dict = {}
    tokeep = set()

    with open("${uk_metadata}", 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

        for row in reader:
            fasta_header = row["fasta_header"]
            id = row["biosample_source_id"]
            epi_day = int(row["edin_epi_day"])

            if id in ["None", "", None]:
                tokeep.add(fasta_header)
                continue

            if id in dup_dict:
                if epi_day < dup_dict[id]["epi_day"]:
                    dup_dict[id] = {"fasta_header": fasta_header, "epi_day": epi_day}
                else:
                    continue
            else:
                dup_dict[id] = {"fasta_header": fasta_header, "epi_day": epi_day}

    for k,v in dup_dict.items():
        tokeep.add(v["fasta_header"])


    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.deduplicated_by_biosamplesourceid.csv", 'w', newline = '') as csv_out, \
         open("${uk_fasta.baseName}.deduplicated_by_biosamplesourceid.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["fasta_header"]
            if fasta_header in tokeep:
                writer.writerow(row)
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\\n")
                fasta_out.write(str(seqrec.seq) + "\\n")
            else:
                continue
    """
}

process uk_remove_duplicates_rootbiosample_by_gaps {
    /**
    * Where duplicate root_biosample, keeps the most complete
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    * @params date
    */

    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.deduplicated_by_rootbiosample.fa", emit: uk_fasta_updated
    path "${uk_metadata.baseName}.deduplicated_by_rootbiosample.csv", emit: uk_metadata_updated

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${uk_fasta}", "fasta")

    dup_dict = {}
    tokeep = set()

    with open("${uk_metadata}", 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

        for row in reader:
            fasta_header = row["fasta_header"]
            id = row["root_biosample_source_id"]
            completeness = float(row["unmapped_genome_completeness"])

            if id in ["None", "", None]:
                tokeep.add(fasta_header)
                continue

            if id in dup_dict:
                if completeness > dup_dict[id]["completeness"]:
                    dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}
                else:
                    continue
            else:
                dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}

    for k,v in dup_dict.items():
        tokeep.add(v["fasta_header"])


    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.deduplicated_by_rootbiosample.csv", 'w', newline = '') as csv_out, \
         open("${uk_fasta.baseName}.deduplicated_by_rootbiosample.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["fasta_header"]
            if fasta_header in tokeep:
                writer.writerow(row)
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\\n")
                fasta_out.write(str(seqrec.seq) + "\\n")
            else:
                continue
    """
}


process unify_headers {
    input:
    path uk_fasta
    path uk_metadata

    output:
    path "${uk_fasta.baseName}.UH.fa"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${uk_fasta}", "fasta")

    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
        open("${uk_fasta.baseName}.UH.fa", "w") as fasta_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        for row in reader:
            record = alignment[row["fasta_header"]]
            fasta_out.write(">" + row["sample_name"] + "\\n")
            fasta_out.write(str(record.seq) + "\\n")
    """
}


workflow deduplicate_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        uk_annotate_with_unmapped_genome_completeness(uk_fasta, uk_metadata)
        uk_remove_duplicates_COGID_by_proportionN(uk_fasta, uk_annotate_with_unmapped_genome_completeness.out)
        uk_remove_duplicates_biosamplesourceid_by_date(uk_remove_duplicates_COGID_by_proportionN.out.uk_fasta_updated, uk_remove_duplicates_COGID_by_proportionN.out.uk_metadata_updated)
        uk_remove_duplicates_rootbiosample_by_gaps(uk_remove_duplicates_biosamplesourceid_by_date.out.uk_fasta_updated, uk_remove_duplicates_biosamplesourceid_by_date.out.uk_metadata_updated)
        uk_unify_headers(uk_remove_duplicates_rootbiosample_by_gaps.out.uk_fasta_updated, uk_remove_duplicates_rootbiosample_by_gaps.out.uk_metadata_updated)
    emit:
        fasta = uk_unify_headers.out
        metadata = uk_remove_duplicates_rootbiosample_by_gaps.out.uk_metadata_updated
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    deduplicate_cog_uk(uk_fasta, uk_metadata)
}
