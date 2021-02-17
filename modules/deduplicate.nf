#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir
publish_dev = file(params.publish_dev)


process annotate_with_unmapped_genome_completeness {
    /**
        * Adds a column to metadata with proportion of genome which is complete
        * @input fasta, metadata
        * @output metadata
        */

        input:
        path fasta
        path metadata

        output:
        path "${metadata.baseName}.annotated.csv"

        script:
        """
        #!/usr/bin/env python3
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index("${fasta}", "fasta")

        with open("${metadata}", 'r', newline = '') as csv_in, \
             open("${metadata.baseName}.annotated.csv", 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["unmapped_genome_completeness"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            id_key = "fasta_header"
            if "edin_header" in reader.fieldnames:
                id_key = "edin_header"

            for row in reader:
                id = row[id_key]
                if id in alignment:
                    seq = str(alignment[id].seq)
                    completeness = float(len(seq.replace("N", "")) / len(seq))
                    row["unmapped_genome_completeness"] = completeness
                    writer.writerow(row)
                else:
                    row["unmapped_genome_completeness"] = 0.0
                    writer.writerow(row)
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
            if row["why_excluded"]:
                continue
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
                if not row["why_excluded"]:
                    row["why_excluded"] = "duplicate central_sample_id"
                writer.writerow(row)
    """
}


process remove_duplicates_by_date {
    /**
    * Where duplicate sequence_name, keeps the earliest
    * @input fasta, metadata
    * @output fasta_updated, metadata_updated
    */

    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.deduplicated.fa", emit: fasta_updated
    path "${metadata.baseName}.deduplicated.csv", emit: metadata_updated

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    dup_dict = {}
    tokeep = set()

    with open("${metadata}", 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

        for row in reader:
            if row["why_excluded"]:
                continue

            fasta_header = row["edin_header"]
            id = row["sequence_name"]
            epi_day = int(row["edin_epi_day"])
            completeness = float(row["unmapped_genome_completeness"])

            if id in ["None", "", None]:
                tokeep.add(fasta_header)
                continue

            if id in dup_dict:
                if epi_day < dup_dict[id][0]["epi_day"]:
                    dup_dict[id].insert(0, {"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness})
                else:
                    dup_dict[id].append({"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness})
            else:
                dup_dict[id] = [{"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness}]

    with open("deduplicated.log", "w") as log:
        for k,v in dup_dict.items():
            tokeep.add(v[0]["fasta_header"])
            if len(v) > 1:
                for dup in v[1:]:
                    log.write("For id %s, %s epi_day:%s completeness:%s kept, %s epi_day:%s completeness:%s removed as duplicate\\n" \
                    %(k, v[0]["fasta_header"], v[0]["epi_day"], v[0]["completeness"], dup["fasta_header"], \
                               dup["epi_day"], dup["completeness"]))

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${metadata}", 'r', newline = '') as csv_in, \
         open("${metadata.baseName}.deduplicated.csv", 'w', newline = '') as csv_out, \
         open("${fasta.baseName}.deduplicated.fa", 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["edin_header"]

            if fasta_header in tokeep:
                writer.writerow(row)
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\\n")
                fasta_out.write(str(seqrec.seq) + "\\n")
            else:
                if not row["why_excluded"]:
                    row["why_excluded"] = "duplicate sequence_name"
                writer.writerow(row)
    """
}


process unify_headers {
    input:
    path fasta
    path metadata

    output:
    path "${fasta.baseName}.UH.fa"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${fasta}", "fasta")

    with open("${metadata}", 'r', newline = '') as csv_in, \
        open("${fasta.baseName}.UH.fa", "w") as fasta_out:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        id_key = "fasta_header"
        if "edin_header" in reader.fieldnames:
            id_key = "edin_header"
        for row in reader:
            if row["why_excluded"]:
                continue
            if row[id_key] in alignment:
                record = alignment[row[id_key]]
                fasta_out.write(">" + row["sequence_name"] + "\\n")
                fasta_out.write(str(record.seq) + "\\n")
    """
}


process uk_label_sourceid_duplicates_to_omit {
    /**
    * Where duplicate source_id, labels all but the earliest as duplicates
    * @input uk_fasta, uk_metadata
    * @output uk_fasta_updated, uk_metadata_updated
    */

    publishDir "${publish_dev}/", pattern: "*.log", mode: 'copy'

    input:
    path uk_metadata

    output:
    path "${uk_metadata.baseName}.deduplicated_by_sourceid.csv", emit: uk_metadata_updated
    path "deduplicated_by_sourceid.log", emit: deduplicate_log

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    dup_dict = {}
    tokeep = set()

    with open("${uk_metadata}", 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

        for row in reader:
            fasta_header = row["sequence_name"]
            id = row["source_id"]
            epi_day = int(row["edin_epi_day"])
            completeness = float(row["unmapped_genome_completeness"])

            if id in ["None", "", None]:
                tokeep.add(fasta_header)
                continue

            if id in dup_dict:
                if epi_day < dup_dict[id][0]["epi_day"]:
                    dup_dict[id].insert(0, {"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness})
                else:
                    dup_dict[id].append({"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness})
            else:
                dup_dict[id] = [{"fasta_header": fasta_header, "epi_day": epi_day, "completeness":completeness}]

    with open("deduplicated_by_sourceid.log", "w") as log:
        for k,v in dup_dict.items():
            tokeep.add(v[0]["fasta_header"])
            if len(v) > 1:
                for dup in v[1:]:
                    log.write("For id %s, %s epi_day:%s completeness:%s kept, %s epi_day:%s completeness:%s removed as duplicate\\n" \
                    %(k, v[0]["fasta_header"], v[0]["epi_day"], v[0]["completeness"], dup["fasta_header"], \
                               dup["epi_day"], dup["completeness"]))


    with open("${uk_metadata}", 'r', newline = '') as csv_in, \
         open("${uk_metadata.baseName}.deduplicated_by_sourceid.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["duplicate"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            row["duplicate"] = None
            fasta_header = row["sequence_name"]
            if fasta_header not in tokeep:
                row["duplicate"] = "True"
            writer.writerow(row)
    """
}


workflow deduplicate_cog_uk {
    take:
        uk_fasta
        uk_metadata
    main:
        annotate_with_unmapped_genome_completeness(uk_fasta, uk_metadata)
        uk_remove_duplicates_COGID_by_proportionN(uk_fasta, annotate_with_unmapped_genome_completeness.out)
        unify_headers(uk_remove_duplicates_COGID_by_proportionN.out.uk_fasta_updated, uk_remove_duplicates_COGID_by_proportionN.out.uk_metadata_updated)
        uk_label_sourceid_duplicates_to_omit(uk_remove_duplicates_COGID_by_proportionN.out.uk_metadata_updated)
    emit:
        fasta = unify_headers.out
        metadata = uk_label_sourceid_duplicates_to_omit.out.uk_metadata_updated
}


workflow deduplicate_gisaid {
    take:
        gisaid_fasta
        gisaid_metadata
    main:
        annotate_with_unmapped_genome_completeness(gisaid_fasta, gisaid_metadata)
        remove_duplicates_by_date(gisaid_fasta, annotate_with_unmapped_genome_completeness.out)
        unify_headers(remove_duplicates_by_date.out.fasta_updated, remove_duplicates_by_date.out.metadata_updated)
    emit:
        fasta = unify_headers.out
        metadata = remove_duplicates_by_date.out.metadata_updated
}


workflow {
    uk_fasta = file(params.uk_fasta)
    uk_metadata = file(params.uk_metadata)
    deduplicate_cog_uk(uk_fasta, uk_metadata)
}
