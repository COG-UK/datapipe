#!/usr/bin/env python3
from Bio import SeqIO
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Split in fasta and metadata into lineageless for pangolin and those with a lineage""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=False, default=None, help='Aligned FASTA')
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV of metadata')
    parser.add_argument('--previous-metadata', dest = 'previous_metadata', required=True, help='CSV of from previous run')
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=False, default=None, help='FASTA to write out')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV of metadata')

    args = parser.parse_args()
    return args

def prepare_for_pangolin(in_fasta, in_metadata, previous_metadata, out_fasta, out_metadata):
    print(in_fasta, in_metadata, previous_metadata, out_fasta, out_metadata)
    if in_fasta:
        alignment = SeqIO.index(in_fasta, "fasta")
    else:
        alignment = None

    taxon = "taxon"
    keys = {"lineage": "lineage",
            "lineages_version": "version",
            "lineage_conflict": "conflict",
            "lineage_ambiguity_score": "ambiguity_score",
            "scorpio_call":"scorpio_call",
            "scorpio_support":"scorpio_support",
            "scorpio_conflict":"scorpio_conflict",
            "usher_lineage":"usher_lineage",
            "usher_lineages_version": "usher_lineages_version"}
    lineage_dict = {}

    with open(previous_metadata, 'r', newline = '') as lineages_in:
        reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")

        if "fasta_header" in reader.fieldnames:
            taxon = "fasta_header"
        elif "edin_header" in reader.fieldnames:
            taxon = "edin_header"
        elif "sequence_name" in reader.fieldnames:
            taxon = "sequence_name"

        if "lineages_version" in reader.fieldnames:
            keys["lineages_version"] = "lineages_version"
        elif "version" in reader.fieldnames:
            keys["lineages_version"] = "version"
        elif "pangoLEARN_version" in reader.fieldnames:
            keys["lineages_version"] =  "pangoLEARN_version"

        for row in reader:
            if row[taxon] in lineage_dict:
                print("%s occurs more than once in lineages input file" % row[taxon])
                continue
            lineage_dict[row[taxon]] = {}
            for key in keys:
                value = keys[key]
                if value in row:
                    lineage_dict[row[taxon]][key] = row[value]


    if out_fasta:
        fasta_out = open(out_fasta, 'w')

    with open(in_metadata, 'r', newline = '') as csv_in, \
        open(out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames
        print(fieldnames, len(fieldnames))
        if len(fieldnames) <= 1:
            csv_in.close()
            csv_in = open(in_metadata, 'r', newline = '')
            reader = csv.DictReader(csv_in, delimiter="\t", quotechar='\"', dialect = "unix")
            fieldnames = reader.fieldnames
        fieldnames.extend([key for key in keys if key not in fieldnames])
        writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        taxon = "taxon"
        if "fasta_header" in reader.fieldnames:
            taxon = "fasta_header"
        elif "edin_header" in reader.fieldnames:
            taxon = "edin_header"
        elif "sequence_name" in reader.fieldnames:
            taxon = "sequence_name"
        print(taxon)
        print(reader.fieldnames)

        missing_lineage = 0

        for row in reader:
            for key in keys:
                if key not in row:
                    row[key] = None

            fasta_header = row[taxon]

            if fasta_header in lineage_dict:
                row.update(lineage_dict[fasta_header])
            elif alignment and fasta_out and fasta_header in alignment:
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\n")
                fasta_out.write(str(seqrec.seq) + "\n")
            if not row["lineage"]:
                missing_lineage += 1
            writer.writerow(row)

    if out_fasta:
        fasta_out.close()

    with open("pango.log", "w") as f:
            f.write("Number of sequences missing lineage assignments: %i" %missing_lineage)

def main():
    args = parse_args()
    print(args)
    prepare_for_pangolin(args.in_fasta, args.in_metadata, args.previous_metadata, args.out_fasta, args.out_metadata)

if __name__ == '__main__':
    main()
