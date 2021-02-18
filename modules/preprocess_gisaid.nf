#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process gisaid_process_json {
    /**
    * Downloads
    * @input json
    * @output gisaid_fasta, gisaid_metadata
    * @params gisaid_omissions
    */

    input:
    path json

    output:
    path "gisaid.fasta", emit: fasta
    path "gisaid.csv", emit: metadata

    script:
    """
    datafunk process_gisaid_data \
        --input-json ${json} \
        --input-metadata False \
        --exclude-file ${gisaid_omissions} \
        --output-fasta "gisaid.fasta" \
        --output-metadata "gisaid.csv" \
        --exclude-undated
    """
}


process gisaid_add_columns_to_metadata {
    input:
    path gisaid_fasta
    path gisaid_metadata

    output:
    path "${gisaid_metadata.baseName}.add_metadata.csv"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import csv

    alignment = SeqIO.index("${gisaid_fasta}", "fasta")

    with open("${gisaid_metadata}", 'r', newline = '') as csv_in, \
        open("${gisaid_metadata.baseName}.add_metadata.csv", 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ['sequence_name', 'why_excluded'], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            edin_header = row["edin_header"]
            new_header = edin_header.split("|")[0]
            row['sequence_name'] = new_header
            if edin_header not in alignment:
                row['why_excluded'] = "filtered during loading from JSON"
            elif row["edin_epi_day"] == '':
                row['why_excluded'] = "no date"
            else:
                row['why_excluded'] = ""
            writer.writerow(row)
    """
}


gisaid_omissions = file(params.gisaid_omissions)

workflow preprocess_gisaid {
    take:
        gisaid_json
    main:
        gisaid_json.splitText( by: params.chunk_size, file: true ).set{ json_chunks }
        gisaid_process_json(json_chunks)
        gisaid_add_columns_to_metadata(gisaid_process_json.out.fasta, gisaid_process_json.out.metadata)
        gisaid_process_json.out.fasta.collectFile(newLine: true).set{ fasta_result }
        gisaid_add_columns_to_metadata.out.collectFile(newLine: true, keepHeader: true, skip: 1)
                                          .set{ metadata_result }
    emit:
        fasta = fasta_result
        metadata = metadata_result
}


workflow {
    gisaid_json = file(params.gisaid_json)

    preprocess_gisaid(gisaid_json)
}
