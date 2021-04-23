#!/usr/bin/env python3

import argparse
import json
import subprocess
import os
import sys
import glob

class Error (Exception): pass

def parse_args():
    parser = argparse.ArgumentParser(description="""Create published files from config file""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--unaligned_fasta', dest = 'unaligned_fasta', required=False, help='Raw FASTA')
    parser.add_argument('--aligned_fasta', dest = 'aligned_fasta', required=False, help='Aligned, masked, untrimmed FASTA')
    parser.add_argument('--trimmed_fasta', dest = 'trimmed_fasta', required=False, help='Aligned, masked, trimmed and filtered FASTA')
    parser.add_argument('--gisaid_fasta', dest = 'global_fasta', required=False, help='GISAID aligned FASTA')
    parser.add_argument('--cog_global_fasta', dest = 'cog_global_fasta', required=False, help='COG GISAID aligned FASTA')

    parser.add_argument('--cog_metadata', dest = 'cog_metadata', required=False, help='MASSIVE CSV')
    parser.add_argument('--gisaid_metadata', dest = 'global_metadata', required=False, help='MASSIVE CSV')
    parser.add_argument('--cog_global_metadata', dest = 'cog_global_metadata', required=False, help='MASSIVE CSV')

    parser.add_argument('--mutations', dest = 'mutations', required=False, help='Mutations CSV')
    parser.add_argument('--constellations', dest = 'constellations', required=False, help='Constellations CSV')

    parser.add_argument('--recipes', dest = 'recipes', required=True, help='JSON of recipes')
    parser.add_argument('--date', dest = 'date', required=True, help='Datestamp for published files')

    args = parser.parse_args()
    return args

#"data": "cog", "gisaid" or "cog_global"
#"fasta": "unaligned", "aligned", "trimmed", "cog_global" or "gisaid"
#"metadata_fields": []
#"mutations": True or False to add columns from mutations
#"shuffle": True to shuffle rows of metadata
#"where": free text to be passed to fastafunk fetch --where-column
#"suffix": something to append to file names
#"exclude_uk": True or False to exclude samples from UK
#"uk_only": True or False to include only samples from UK from cog_global metadata
#"drop_index": name of index column that should be dropped at the end

def get_info_from_config(config_dict, outdir, date, fasta_dict, csv_dict, mutations_file, constellations_file):
    info_dict = {"suffix":None, "data":None, "fasta":None, "metadata_fields":None,
                 "where": None, "mutations":False, "constellations":False, 
                 "shuffle":False, "drop_index": None,
                 "exclude_uk":False, "uk_only": False, "exclude_cog":False, "cog_only": False,
                 "date": date,
                 "in_fa":None, "in_csv":None, "in_muts":None, "in_con":None,
                 "out_fa":"tmp.fa", "intermediate_csv":"tmp.csv", "out_csv":"tmp.csv"}
    info_dict.update(config_dict)

    if info_dict["fasta"] in fasta_dict.keys():
        info_dict["in_fa"] = fasta_dict[info_dict["fasta"]]
    elif info_dict["data"] == "cog_global":
        info_dict["in_fa"] = fasta_dict["cog_global"]
    elif info_dict["data"] == "gisaid":
        info_dict["in_fa"] = fasta_dict["gisaid"]
    elif info_dict["data"] == "cog":
            info_dict["in_fa"] = fasta_dict["trimmed"]
    else:
        sys.exit("Config entries need to specify either fasta in ['unaligned', 'aligned', 'trimmed', 'cog_global', 'gisaid'] or data \
        in ['cog', 'cog_global', 'gisaid']")

    if info_dict["data"] is None:
        if info_dict["fasta"] == "cog_global":
            info_dict["data"] = "cog_global"
        elif info_dict["fasta"] == "gisaid":
            info_dict["data"] = "gisaid"
        else:
            info_dict["data"] = "cog"

    if info_dict["data"] == "cog_global":
            info_dict["in_csv"] = csv_dict["cog_global"]
    elif info_dict["data"] == "cog":
            info_dict["in_csv"] = csv_dict["cog"]
    elif info_dict["data"] == "gisaid":
            info_dict["in_csv"] = csv_dict["gisaid"]

    info_dict["in_muts"] = mutations_file
    info_dict["in_con"] = constellations_file

    start = "%s/%s_%s" %(outdir, info_dict["data"], info_dict["date"])
    if info_dict["suffix"]:
        start += "_%s" %info_dict["suffix"]
    csv_end = ".csv"

    if info_dict["fasta"]:
        csv_end = "_metadata.csv"
        if info_dict["fasta"]=="aligned" or (info_dict["metadata_fields"] and info_dict["fasta"]!="unaligned"):
            info_dict["out_fa"] = "%s_alignment.fa" %start
        else:
            info_dict["out_fa"] = "%s.fa" %start

    info_dict["out_csv"] = "%s%s" %(start, csv_end)

    if info_dict["out_fa"] != "tmp.fa" and info_dict["in_fa"] is None:
        sys.exit("Please provide the appropriate FASTA file")
    if info_dict["metadata_fields"] is not None and info_dict["in_csv"] is None:
        sys.exit("Please provide the appropriate CSV file")
    if info_dict["mutations"] is not None and info_dict["in_muts"] is None:
        sys.exit("Please provide the appropriate mutations file")
    if info_dict["constellations"] is not None and info_dict["in_con"] is None:
            sys.exit("Please provide the appropriate mutations file")
    print(info_dict)
    return info_dict

def syscall(cmd_list, allow_fail=False):
    if None in cmd_list:
        print('None in list', cmd_list, file=sys.stderr)
        raise Error('Error in command. Cannot continue')
    command = ' '.join(cmd_list)
    print(command)
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if (not allow_fail) and completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('\nOutput from stdout:', completed_process.stdout, sep='\n', file=sys.stderr)
        print('\nOutput from stderr:', completed_process.stderr, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')
    print(completed_process.stdout)
    return completed_process

def publish_file(outdir, info_dict):
    if info_dict["metadata_fields"] is None:
        cmd_list = ["cp", info_dict["in_fa"], info_dict["out_fa"]]
        syscall(cmd_list)
        return

    if info_dict["exclude_uk"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.no_uk.csv --column is_uk --is_true"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.no_uk.csv"

    if info_dict["exclude_cog"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.no_cog.csv --column is_cog_uk --is_true"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.no_cog.csv"

    if info_dict["uk_only"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.uk_only.csv --column is_uk --is_false"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.uk_only.csv"

    if info_dict["cog_only"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.cog_only.csv --column is_cog_uk --is_false"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.cog_only.csv"

    if info_dict["shuffle"]:
        cmd_list = ["fastafunk shuffle --in-metadata", info_dict["in_csv"], "--out-metadata", "tmp.shuffled.csv"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.shuffled.csv"

    cmd_list = ["fastafunk fetch --in-fasta", info_dict["in_fa"], "--in-metadata", info_dict["in_csv"],
              "--index-column sequence_name --out-fasta", info_dict["out_fa"],
              "--out-metadata", info_dict["intermediate_csv"], "--restrict --low-memory --keep-omit-rows"]
    if info_dict["metadata_fields"]:
            cmd_list.append("--filter-column")
            cmd_list.extend(info_dict["metadata_fields"])
    if info_dict["where"]:
        cmd_list.append("--where-column %s" %info_dict["where"])
    syscall(cmd_list)

    if info_dict["mutations"]:
        cmd_list = ["fastafunk add_columns --in-metadata", info_dict["intermediate_csv"],
        "--in-data", info_dict["in_muts"], "--index-column sequence_name",
        "--join-on sequence_name --out-metadata tmp.muts.csv"]
        info_dict["intermediate_csv"] = "tmp.muts.csv"
        syscall(cmd_list)

    if info_dict["constellations"]:
            cmd_list = ["fastafunk add_columns --in-metadata", info_dict["intermediate_csv"],
            "--in-data", info_dict["in_con"], "--index-column sequence_name",
            "--join-on sequence_name --out-metadata tmp.constellations.csv"]
            info_dict["intermediate_csv"] = "tmp.constellations.csv"
            syscall(cmd_list)

    if info_dict["drop_index"]:
        cmd_list = ["fastafunk drop_columns --in-metadata", info_dict["intermediate_csv"],
        "--columns", info_dict["drop_index"],
        "--out-metadata tmp.anon.csv"]
        info_dict["intermediate_csv"] = "tmp.anon.csv"
        syscall(cmd_list)


    cmd_list = ["mv", info_dict["intermediate_csv"], info_dict["out_csv"]]
    syscall(cmd_list)

    #tmp = glob.glob("tmp.*")
    #if len(tmp) > 0:
    #    cmd_list = ["rm tmp.*"]
    #    syscall(cmd_list)

def main():
    args = parse_args()
    print(args)
    fasta_dict = {"unaligned":args.unaligned_fasta, "aligned":args.aligned_fasta, "trimmed":args.trimmed_fasta, "cog_global": args.cog_global_fasta, "gisaid": args.global_fasta}
    print(fasta_dict)
    csv_dict = {"cog":args.cog_metadata, "cog_global":args.cog_global_metadata, "gisaid": args.global_metadata}
    print(csv_dict)
    mutations_file = args.mutations
    print(mutations_file)
    constellations_file = args.constellations
    print(constellations_file)

    recipes = {}
    with open(args.recipes, 'r') as f:
        recipes = json.load(f)

    for outdir in recipes.keys():
        os.makedirs(outdir,exist_ok=True)
        for recipe in recipes[outdir]:
            info_dict = get_info_from_config(recipe, outdir, args.date, fasta_dict, csv_dict, mutations_file, constellations_file)
            publish_file(outdir, info_dict)

if __name__ == '__main__':
    main()
