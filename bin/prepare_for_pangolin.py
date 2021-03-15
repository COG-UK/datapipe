#!/usr/bin/env python3
from Bio import SeqIO
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Split in fasta and metadata into lineageless for pangolin and those with a lineage""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Aligned FASTA')
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV of metadata')
    parser.add_argument('--previous-metadata', dest = 'previous_metadata', required=True, help='CSV of from previous run')
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=True, help='FASTA to write out')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV of metadata')

    args = parser.parse_args()
    return args

def prepare_for_pangolin(in_fasta, in_metadata, previous_metadata, out_fasta, out_metadata):
    alignment = SeqIO.index(in_fasta, "fasta")

    lineage_dict = {}
    (taxon,lin,pango_version,prob) = "taxon","lineage","pangoLEARN_version","probability"
    with open(previous_metadata, 'r', newline = '') as lineages_in:
        reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")

        if "fasta_header" in reader.fieldnames:
            taxon = "fasta_header"
        elif "edin_header" in reader.fieldnames:
            taxon = "edin_header"
        elif "sequence_name" in reader.fieldnames:
            taxon = "sequence_name"
        if "lineages_version" in reader.fieldnames:
            pango_version = "lineages_version"
        if "lineage_support" in reader.fieldnames:
            prob = "lineage_support"

        for row in reader:
            if row[taxon] in lineage_dict:
                print("%s occurs more than once in lineages input file" % row[taxon])
                continue
            lineage_dict[row[taxon]] = {"lineage": row[lin], "pangoLEARN_version": row[pango_version], "probability": row[prob]}

    with open(in_metadata, 'r', newline = '') as csv_in, \
        open(out_metadata, 'w', newline = '') as csv_out, \
        open(out_fasta, 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["lineage", "pangoLEARN_version", "probability"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        if "fasta_header" in reader.fieldnames:
            taxon = "fasta_header"
        elif "edin_header" in reader.fieldnames:
            taxon = "edin_header"

        for row in reader:
            lineage = ""
            pangoLEARN_version = ""
            probability = ""

            fasta_header = row[taxon]

            if fasta_header in lineage_dict:
                lineage = lineage_dict[fasta_header]["lineage"]
                pangoLEARN_version = lineage_dict[fasta_header]["pangoLEARN_version"]
                probability = lineage_dict[fasta_header]["probability"]
            elif fasta_header in alignment:
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\n")
                fasta_out.write(str(seqrec.seq) + "\n")

            row["lineage"] = lineage
            row["pangoLEARN_version"] = pangoLEARN_version
            row["probability"] = probability

            writer.writerow(row)

def main():
    args = parse_args()
    prepare_for_pangolin(args.in_fasta, args.in_metadata, args.previous_metadata, args.out_fasta, args.out_metadata)

if __name__ == '__main__':
    main()
