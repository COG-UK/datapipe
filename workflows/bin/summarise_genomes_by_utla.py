#!/usr/bin/env python3
import csv
from collections import defaultdict
from collections import Counter
import datetime as dt

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", action="store")
parser.add_argument("--date", action='store')
args = parser.parse_args()

metadata = args.metadata
file_date = args.date

def main(metadata, file_date):

    utla_to_region = {}
    utla_to_code = {}

    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['utla'] != "" and "|" not in l['utla']:
                utla_to_region[l['utla']] = l['NUTS1']
                utla_to_code[l['utla']] = l['utla_code']

    utla_delta = defaultdict(list)
    utla_other = defaultdict(list)
    utla_all = defaultdict(list)
    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['sample_date'] != "":
                date = dt.datetime.strptime(l['sample_date'],"%Y-%m-%d").date()
                if l['utla'] != "" and "|" not in l['utla']:
                    if l['scorpio_call'] == "Delta (B.1.617.2-like)":
                        utla_delta[date].append(l['utla'])
                    else:
                        utla_other[date].append(l['utla'])

                    utla_all[date].append(l['utla'])

    delta_counts = {}
    other_counts = {}
    all_counts = {}

    for k,v in utla_delta.items():
        delta_counts[k] = Counter(v)

    for k,v in utla_other.items():
        other_counts[k] = Counter(v)

    for k,v in utla_all.items():
        all_counts[k] = Counter(v)

    fieldnames = ["date", "utla", "utla_code", "NUTS1", "delta_count", "other_count", "total_count"]
    with open(f"UTLA_genome_counts_{file_date}.csv", 'w') as fw:
        writer = csv.DictWriter(fw, fieldnames=fieldnames)
        writer.writeheader()
        for date, utla_dict in sorted(all_counts.items()):
            for utla, count in utla_dict.items():
                write_dict = {}
                write_dict["date"] = date
                write_dict["utla"] = utla
                write_dict["utla_code"] = utla_to_code[utla]
                write_dict["NUTS1"] = utla_to_region[utla]
                write_dict["total_count"] = count
                if date in delta_counts:
                    if utla in delta_counts[date]:
                        write_dict["delta_count"] = delta_counts[date][utla]
                    else:
                        write_dict["delta_count"] = 0
                else:
                    write_dict["delta_count"] = 0

                if date in other_counts:
                    if utla in other_counts[date]:
                        write_dict["other_count"] = other_counts[date][utla]
                    else:
                        write_dict["other_count"] = 0
                else:
                    write_dict["other_count"] = 0
                    
                writer.writerow(write_dict)


if __name__ == '__main__':
    main(metadata, file_date)