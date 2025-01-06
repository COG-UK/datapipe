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
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=True, help='FASTA to write out')
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='FASTA')

    args = parser.parse_args()

    return args

def run(in_fasta, in_metadata, out_fasta, out_metadata):
    dup_dict = {}
    tokeep = set()

    with open(in_metadata, 'r', newline = '') as csv_in:
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
                    log.write("For id %s, %s epi_day:%s completeness:%s kept, %s epi_day:%s completeness:%s removed as duplicate\n" \
                    %(k, v[0]["fasta_header"], v[0]["epi_day"], v[0]["completeness"], dup["fasta_header"], \
                               dup["epi_day"], dup["completeness"]))

    alignment = SeqIO.index(in_fasta, "fasta")

    with open(in_metadata, 'r', newline = '') as csv_in, \
         open(out_metadata, 'w', newline = '') as csv_out, \
         open(out_fasta, 'w') as fasta_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["edin_header"]

            if fasta_header in tokeep:
                writer.writerow(row)
                seqrec = alignment[fasta_header]
                fasta_out.write(">" + seqrec.id + "\n")
                fasta_out.write(str(seqrec.seq) + "\n")
            else:
                if not row["why_excluded"]:
                    row["why_excluded"] = "duplicate sequence_name"
                writer.writerow(row)

def main():
    args = parse_args()
    run(args.in_fasta, args.in_metadata, args.out_fasta, args.out_metadata)

if __name__ == '__main__':
    main()
