#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Add sample_date, pillar_2 and sequence_name columns""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='TSV from MAJORA')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--accession-file', dest = 'accession_file', required=True, help='TSV of accessions')
    parser.add_argument('--log-file', dest = 'log_file', required=False, help='Log file')

    args = parser.parse_args()

    return args

def add_sample_date(df):
    sample_date = []

    for i,row in df.iterrows():
        if not pd.isnull(row['collection_date']) and row['collection_date'] != "None":
            sample_date.append(row['collection_date'])
        elif not pd.isnull(row['received_date']) and row['received_date'] != "None":
            sample_date.append(row['received_date'])
        else:
            sample_date.append("")

    df['sample_date'] = sample_date
    return df

def add_pillar_2(df):
    pillar_2 = []

    for i,row in df.iterrows():
        if row['collection_pillar'] == 2 or row['central_sample_id'][0:4] in ["ALDP", "CAMC", "MILK", "QEUH"]:
            pillar_2.append(True)
        else:
            pillar_2.append(False)

    df['pillar_2'] = pillar_2
    return df

def add_sequence_name(df):
    adm1_to_country = {"UK-SCT": "Scotland",
                       "UK-WLS": "Wales",
                       "UK-ENG": "England",
                       "UK-NIR": "Northern_Ireland"}

    sequence_name = []

    for i,row in df.iterrows():
        country = adm1_to_country[row['adm1']]
        id = row['central_sample_id']
        year = str(row['sample_date']).split("-")[0]
        name = country + "/" + id + "/" + year

        sequence_name.append(name)

    df['sequence_name'] = sequence_name
    return df

def load_accessions(accession_file, log_handle):
    df_acc = pd.read_csv(accession_file, sep='\t')
    accessions_dict = {}
    for i,row in df_acc.iterrows():
        central_sample_id = row["central_sample_id"]
        run_name = row["run_name"]
        gisaid_accession = row["gisaid.accession"]

        if central_sample_id in accessions_dict:
            if run_name in accessions_dict[central_sample_id]:
                log_handle.write(f'duplicate central_sample_id * run_name in accessions list: {central_sample_id} {run_name}\n')
                continue
            accessions_dict[central_sample_id][run_name] = gisaid_accession
        else:
            accessions_dict[central_sample_id] = {run_name: gisaid_accession}
    return accessions_dict

def add_covv_accession_id(df, accession_file, log_handle):
    accessions_dict = load_accessions(accession_file, log_handle)

    covv_accession_id = []
    for i,row in df.iterrows():
        acc = ""
        if row["central_sample_id"] in accessions_dict:
            if row["run_name"] in accessions_dict[row["central_sample_id"]]:
                acc = accessions_dict[row["central_sample_id"]][row["run_name"]]

        covv_accession_id.append(acc)

    df['covv_accession_id'] = covv_accession_id
    return df

def main():
    args = parse_args()
    if args.log_file:
        log_handle = open(args.log_file, 'w')
    else:
        log_handle = sys.stdout

    df = pd.read_csv(args.in_metadata, sep = "\t")
    df = add_sample_date(df)
    df = add_pillar_2(df)
    df = add_sequence_name(df)
    df = add_covv_accession_id(df, args.accession_file, log_handle)

    df.to_csv(args.out_metadata, index=False, sep = ",")

    log_handle.close()

if __name__ == '__main__':
    main()
