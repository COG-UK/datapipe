#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

project_dir = projectDir


process uk_strip_header_digits {
    /**
    * Strips extra header info from FASTA
    * @input uk_fasta
    * @output uk_fasta_updated
    */

    input:
    file uk_fasta

    output:
    file "${uk_fasta.baseName}.updated.fasta"

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    fasta_in = SeqIO.parse("${uk_fasta}", "fasta")
    with open("${uk_fasta.baseName}.updated.fasta", 'w') as f:
        for record in fasta_in:
            ID = record.description.split("|")[0]
            f.write(">" + ID + "\\n")
            f.write(str(record.seq) + "\\n")
    """
}


process uk_add_columns_to_metadata {
    /**
    * Takes the MAJORA TSV of metadata and adds/updates columns for sample_date, pillar_2,
    * sequence_name and covv_accession_id
    * @input uk_metadata, uk_accessions
    * @output uk_metadata_updated
    */

    input:
    file uk_metadata
    file uk_accessions

    output:
    file "${uk_metadata.baseName}.updated.csv"

    script:
    """
    $project_dir/../bin/add_to_uk_metadata.py \
        --in-metadata ${uk_metadata} \
        --out-metadata ${uk_metadata.baseName}.updated.csv \
        --accession-file ${uk_accessions}
    """
}


process uk_annotate_to_remove_duplicates {

    input:
    file uk_fasta
    file uk_metadata

    output:
    file "${uk_metadata.baseName}.annotated.csv"

    script:
    """
    fastafunk annotate \
        --in-fasta ${uk_fasta} \
        --in-metadata ${uk_metadata} \
        --out-metadata ${uk_metadata.baseName}.annotated.csv \
        --index-column fasta_header
    """
}


workflow reformat_fasta_and_metadata {
    take:
        uk_fasta
        uk_metadata
        uk_accessions
    main:
        uk_strip_header_digits(uk_fasta)
        uk_add_columns_to_metadata(uk_metadata, uk_accessions)
        uk_annotate_to_remove_duplicates(uk_strip_header_digits.out, uk_add_columns_to_metadata.out)
    emit:
        uk_strip_header_digits.out
        uk_annotate_to_remove_duplicates.out
}

params.uk_fasta = file("test/matched.fa")
params.uk_metadata = file("test/matched.tsv")
params.uk_accessions = file("test/accessions.tsv")

workflow {
    reformat_fasta_and_metadata(params.uk_fasta, params.uk_metadata, params.uk_accessions)
}
