#!/usr/bin/env python3

import argparse
import json
import subprocess
import os
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="""Create published files from config file""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--unaligned_fasta', dest = 'unaligned_fasta', required=True, help='Raw FASTA')
    parser.add_argument('--aligned_fasta', dest = 'aligned_fasta', required=True, help='Aligned, masked, untrimmed FASTA')
    parser.add_argument('--trimmed_fasta', dest = 'trimmed_fasta', required=True, help='Aligned, masked, trimmed and filtered FASTA')
    parser.add_argument('--cog_global_fasta', dest = 'cog_global_fasta', required=True, help='COG GISAID aligned FASTA')

    parser.add_argument('--cog_metadata', dest = 'cog_metadata', required=True, help='MASSIVE CSV')
    parser.add_argument('--cog_global_metadata', dest = 'cog_global_metadata', required=True, help='MASSIVE CSV')

    parser.add_argument('--cog_variants', dest = 'cog_variants', required=True, help='Mutations CSV')
    parser.add_argument('--cog_global_variants', dest = 'cog_global_variants', required=True, help='Mutations CSV')

    parser.add_argument('--recipes', dest = 'recipes', required=True, help='JSON of recipes')
    parser.add_argument('--date', dest = 'date', required=True, help='Datestamp for published files')

    args = parser.parse_args()
    return args

#"data": "cog" or "cog_global"
#"fasta": "unaligned", "aligned", "trimmed", "cog_global"
#"metadata_fields": []
#"variants": True or False to add columns from variants
#"where": free text to be passed to fastafunk fetch --where-column
#"suffix": something to append to file names

def get_info_from_config(config_dict, outdir, date, fasta_dict, csv_dict, var_dict):
    info_dict = {"suffix":None, "data":None, "fasta":None, "metadata_fields":None,
                 "where": None, "variants":False, "date": date,
                 "in_fa":None, "in_csv":None, "in_var":None,
                 "out_fa":"tmp.fa", "out_csv":"tmp.csv", "out_var":None}
    info_dict.update(config_dict)

    if info_dict["fasta"] in fasta_dict.keys():
        info_dict["in_fa"] = fasta_dict[info_dict["fasta"]]
    elif info_dict["data"] == "cog_global":
        info_dict["in_fa"] = fasta_dict["cog_global"]
    elif info_dict["data"] == "cog":
            info_dict["in_fa"] = fasta_dict["trimmed"]
    else:
        sys.exit("Config entries need to specify either fasta in ['unaligned', 'aligned', 'trimmed', 'cog_global'] or data \
        in ['cog', 'cog_global']")

    if info_dict["data"] == "cog_global":
            info_dict["in_csv"] = csv_dict["cog_global"]
            info_dict["in_var"] = var_dict["cog_global"]
    elif info_dict["data"] == "cog":
            info_dict["in_csv"] = csv_dict["cog"]
            info_dict["in_var"] = var_dict["cog"]

    if info_dict["data"] is None:
        if info_dict["fasta"] == "cog_global":
            info_dict["data"] = "cog_global"
        else:
            info_dict["data"] = "cog"

    if info_dict["variants"]:
        if info_dict["suffix"] is None:
            info_dict["out_var"] = "%s/%s_%s_metadata.csv" %(outdir, info_dict["data"], info_dict["date"])
        else:
            info_dict["out_var"] = "%s/%s_%s_%s_metadata.csv" %(outdir, info_dict["data"], info_dict["date"], info_dict["suffix"])
    elif info_dict["metadata_fields"]:
        if info_dict["suffix"] is None:
            info_dict["out_csv"] = "%s/%s_%s_metadata.csv" %(outdir, info_dict["data"], info_dict["date"])
        else:
            info_dict["out_csv"] = "%s/%s_%s_%s_metadata.csv" %(outdir, info_dict["data"], info_dict["date"], info_dict["suffix"])

    if info_dict["fasta"]:
        if info_dict["metadata_fields"]:
            if info_dict["suffix"] is None:
                info_dict["out_fa"] = "%s/%s_%s_alignment.fa" %(outdir, info_dict["data"], info_dict["date"])
            else:
                info_dict["out_fa"] = "%s/%s_%s_%s_alignment.fa" %(outdir, info_dict["data"], info_dict["date"], info_dict["suffix"])
        else:
            if info_dict["suffix"] is None:
                info_dict["out_fa"] = "%s/%s_%s.fa" %(outdir, info_dict["data"], info_dict["date"])
            else:
                info_dict["out_fa"] = "%s/%s_%s_%s.fa" %(outdir, info_dict["data"], info_dict["date"], info_dict["suffix"])

    return info_dict


def publish_file(outdir, info_dict):
    if info_dict["metadata_fields"] is None:
        cmd_list = ["cp", info_dict["in_fa"], info_dict["out_fa"]]
        subprocess.run(' '.join(cmd_list), shell=True)
        return

    cmd_list = ["fastafunk fetch --in-fasta", info_dict["in_fa"], "--in-metadata", info_dict["in_csv"],
              "--index-column sequence_name --out-fasta", info_dict["out_fa"],
              "--out-metadata", info_dict["out_csv"], "--restrict"]
    if info_dict["where"]:
        cmd_list.append("--where-column %s" %info_dict["where"])
    subprocess.run(' '.join(cmd_list), shell=True)

    if info_dict["variants"]:
        cmd_list = ["fastafunk add_columns --in-metadata", info_dict["out_csv"],
        "--in-data", info_dict["in_var"], "--index-column sequence_name",
        "--join-on query --out-metadata", info_dict["out_var"]]
        subprocess.run(' '.join(cmd_list), shell=True)

    tmp = glob.glob("tmp.*")
    if len(tmp) > 0:
        subprocess.run("rm tmp.*", shell=True)

def main():
    args = parse_args()

    fasta_dict = {"unaligned":args.unaligned_fasta, "aligned":args.aligned_fasta, "trimmed":args.trimmed_fasta, "cog_global": args.cog_global_fasta}
    csv_dict = {"cog":args.cog_metadata, "cog_global":args.cog_global_metadata}
    var_dict = {"cog":args.cog_variants, "cog_global":args.cog_global_variants}

    recipes = {}
    with open(args.recipes, 'r') as f:
        recipes = json.load(f)

    for outdir in recipes.keys():
        os.makedirs(outdir,exist_ok=True)
        for recipe in recipes[outdir]:
            info_dict = get_info_from_config(recipe, outdir, args.date, fasta_dict, csv_dict, var_dict)
            publish_file(outdir, info_dict)

if __name__ == '__main__':
    main()
