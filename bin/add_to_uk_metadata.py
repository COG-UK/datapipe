#!/usr/bin/env python3

import sys
import argparse
import csv
from itertools import chain
from epiweeks import Week,Year
from datetime import datetime

adm1a_to_country = {"UK-SCT": "Scotland",
                    "UK-WLS": "Wales",
                    "UK-ENG": "England",
                    "UK-NIR": "Northern_Ireland"}

def parse_args():
    parser = argparse.ArgumentParser(description="""Add sample_date, pillar_2 and sequence_name columns""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='TSV from MAJORA')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--accession-file', dest = 'accession_file', required=False, help='TSV of accession')
    parser.add_argument('--updated-date-file', dest = 'updated_date_file', required=False, help='CSV of date corrections')
    parser.add_argument('--log-file', dest = 'log_file', required=False, help='Log file')

    args = parser.parse_args()

    return args

def load_updated_dates(updated_date_file):
    date_dict = {}
    if updated_date_file:
        with open(updated_date_file, 'r', newline = '') as dates_in:
            reader = csv.DictReader(dates_in, delimiter=",", quotechar='\"', dialect = "unix")
            for row in reader:
                date_dict[row["central_sample_id"]] = row["sample_date"]
    return date_dict

def add_sample_date(row, date_dict):
    if row["central_sample_id"] in date_dict:
        row["sample_date"] = date_dict[row["central_sample_id"]]
        return
    try:
        date = datetime.strptime(row["collection_date"], '%Y-%m-%d').date()
        row["sample_date"] = row["collection_date"]
    except:
        try:
            date = datetime.strptime(row["received_date"], '%Y-%m-%d').date()
            row["sample_date"] = row["received_date"]
        except:
            row["sample_date"] = ""

def add_pillar_2(row):
    if row['collection_pillar'] == 2 or row['central_sample_id'][0:4] in ["ALDP", "CAMC", "MILK", "QEUH"]:
        row["pillar_2"] = True
    else:
        row["pillar_2"] = False

def add_sequence_name(row):
    country = adm1a_to_country[row['adm1']]
    id = row['central_sample_id']
    year = str(row['sample_date']).split("-")[0]
    name = country + "/" + id + "/" + year

    row["sequence_name"] = name

def load_accession(accession_file, log_handle):
    if not accession_file:
        return {}

    accession_dict = {}

    with open(str(accession_file), 'r', newline = '') as acc_in:
        reader = csv.DictReader(acc_in, delimiter="\t", quotechar='\"', dialect = "unix")
        for row in reader:
            central_sample_id = row["central_sample_id"]
            run_name = row["run_name"]
            gisaid_accession = row["gisaid.accession"]

            if central_sample_id in accession_dict:
                if run_name in accession_dict[central_sample_id]:
                    log_handle.write(f'duplicate central_sample_id * run_name in accession list: {central_sample_id} {run_name}\n')
                    continue
                accession_dict[central_sample_id][run_name] = gisaid_accession
            else:
                accession_dict[central_sample_id] = {run_name: gisaid_accession}
    return accession_dict

def add_covv_accession_id(row, accession_dict):
    acc = ""
    if row["central_sample_id"] in accession_dict:
        if row["run_name"] in accession_dict[row["central_sample_id"]]:
            acc = accession_dict[row["central_sample_id"]][row["run_name"]]

    row["covv_accession_id"] = acc

def date_string_to_epi_week(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi week which is cumulative total epidemiological
    weeks since 2019-12-22. Week beginning 2019-12-22 is week 0
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week:
    week = Week.fromdate(date)
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    elif week.year == 2019:
        return("0")
    else:
        cum_epi_week = week.week + len(list(chain(*[[x for x in Year(y).iterweeks()] for y in range(2020, week.year)])))
        return str(cum_epi_week)

def date_string_to_epi_day(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi day which is cumulative total days since 2019-12-22
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week week:
    week = Week.fromdate(date)
    # this is day 1 of epi-week 0:
    day_one = datetime.strptime("2019-12-22", '%Y-%m-%d').date()
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    else:
        cum_epi_day = (date - day_one).days + 1
        return str(cum_epi_day)

def add_epi_week_and_day(row):
    date_str = row["sample_date"]
    epi_week = date_string_to_epi_week(date_str)
    epi_day = date_string_to_epi_day(date_str)

    row["edin_epi_week"] = epi_week
    row["edin_epi_day"] = epi_day

def United_Kingdom_to_UK(row):
    row["adm0"] = row["adm0"].replace("United Kingdom", "UK")

def main():
    args = parse_args()
    if args.log_file:
        log_handle = open(args.log_file, 'w')
    else:
        log_handle = sys.stdout

    date_dict = load_updated_dates(args.updated_date_file)
    accession_dict = load_accession(args.accession_file, log_handle)
    new_columns = ["sample_date", "pillar_2", "sequence_name", "covv_accession_id", "edin_epi_week", "edin_epi_day", "why_excluded"]

    with open(args.in_metadata, 'r', newline = '') as csv_in, \
         open(args.out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter="\t", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + new_columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            try:
                add_sample_date(row, date_dict)
                add_pillar_2(row)
                add_sequence_name(row)
                add_covv_accession_id(row, accession_dict)
                add_epi_week_and_day(row)
                United_Kingdom_to_UK(row)
                row["why_excluded"] = ""
                writer.writerow(row)
            except:
                log_handle.write(f"Error updating metadata for row")
                log_handle.write(str(row))
                sys.exit("Could not update metadata for row, check metadata fields")


    log_handle.close()

if __name__ == '__main__':
    main()
