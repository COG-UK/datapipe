#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="""Add sample_date, is_pillar_2 and sequence_name columns""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='TSV from MAJORA')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='FASTA')

    args = parser.parse_args()

    return args

def run(in_fasta, in_metadata, out_metadata):
    alignment = SeqIO.index(in_fasta, "fasta")

    with open(in_metadata, 'r', newline = '') as csv_in, \
        open(out_metadata, 'w', newline = '') as csv_out:

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
                if len(seq) == 0:
                    print(id)
                    row["unmapped_genome_completeness"] = 0.0
                else:
                    completeness = float(len(seq.replace("N", "")) / len(seq))
                    row["unmapped_genome_completeness"] = completeness
                writer.writerow(row)
            else:
                row["unmapped_genome_completeness"] = 0.0
                writer.writerow(row)

def main():
    args = parse_args()
    run(args.in_fasta, args.in_metadata, args.out_metadata)

if __name__ == '__main__':
    main()
