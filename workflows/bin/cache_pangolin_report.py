#!/usr/bin/env python3

import argparse
import csv
from Bio import SeqIO
import hashlib


def parse_args():
    parser = argparse.ArgumentParser(description="""Add sample_date, is_pillar_2 and sequence_name columns""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='Lineage report from pangolin')
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Unaligned fasta')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='Hashed lineage report from pangolin')

    args = parser.parse_args()

    return args

def get_hash_string(record):
    seq = str(record.seq).upper().encode()
    hash_object = hashlib.md5(seq)
    hash_string = hash_object.hexdigest()
    return hash_string

def cache_report(in_fasta, in_metadata, out_metadata):
    hashed_seqs = set()
    records = SeqIO.index(in_fasta, "fasta")
    index_column = "taxon"

    with open(in_metadata, 'r', newline = '') as csv_in, \
         open(out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames[:]
        fieldnames.remove(index_column)
        print(fieldnames)
        writer = csv.DictWriter(csv_out, fieldnames = ["hash"] + fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            print(row)
            if row[index_column] not in records:
                continue
            record = records[row[index_column]]
            hash = get_hash_string(record)
            if hash not in hashed_seqs:
                del row[index_column]
                row["hash"] = hash
                hashed_seqs.add(hash)
                writer.writerow(row)


def main():
    args = parse_args()
    cache_report(args.in_fasta, args.in_metadata, args.out_metadata)


if __name__ == '__main__':
    main()
