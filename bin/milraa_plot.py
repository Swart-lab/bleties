#!/usr/bin/env python3

import re
import json
import argparse
import sys
import os

import pandas as pd
import matplotlib as mpl
mpl.use("Agg") # allow run without X-server
import matplotlib.pyplot as plt

from bleties import SharedFunctions


# Argument parser
parser = argparse.ArgumentParser(
        description="Plot IES statistics from MILRAA output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gff",
        help="GFF file produced by BleTIES MILRAA module")
parser.add_argument("--type", default='subreads',
        help="Type of reads used to map, either 'ccs' or 'subreads'")
parser.add_argument("--out", "-o", default="test",
        help="""
        Prefix for output files. Three histograms are produced with the
        following filename suffixes: 
            - ies_retention_score
            - ies_length_distribution
            - ies_length_distribution_detail
        """)
parser.add_argument("--out_fmt", default="png",
        help="Format for output files, passed to matplotlib: 'png','pdf','svg'")
parser.add_argument("--hist_len_min", default=26, type=int,
        help="Minimum length to show in histogram of IES lengths (detail)")
parser.add_argument("--hist_len_max", default=400, type=int,
        help="Maximum length to show in histogram of IES lengths (detail)")
parser.add_argument("--hist_style", default="facet",
        help="Style for histograms of IES lengths: 'facet', 'bar', or 'barstacked'")
args = parser.parse_args()

# Import GFF records
gff_obj = SharedFunctions.Gff()
gff_obj.file2gff(args.gff)

# Get relevant fields
out = []
for gff_id in gff_obj:
    row = {}
    row['id'] = gff_id
    if int(gff_obj.getValue(gff_id, "start")) == int(gff_obj.getValue(gff_id, "end")):
        row['type'] = "ins"
    elif int(gff_obj.getValue(gff_id, "start")) < int(gff_obj.getValue(gff_id, "end")):
        row['type'] = "del"
    else:
        raise Exception(f"Invalid GFF3, start > end: check {row['id']}")
    # parse differently for subreads vs. ccs
    if args.type == 'ccs':
        row['mean_cov'] = int(gff_obj.getAttr(gff_id, "average_coverage"))
        # TODO: remove this cigar total? used for plotting and is experimental
        cigar = gff_obj.getAttr(gff_id, "cigar")
        cigar_covs = re.findall(r"\d+[ID]\*(\d+)", cigar)
        cigar_total = sum([int(i) for i in cigar_covs])
        if row['type'] == "ins":
            iesplus = cigar_total
        elif row['type'] == "del":
            iesplus = row['mean_cov'] - cigar_total
        row['retention_score'] = round(float(iesplus / row['mean_cov']), 4)
    elif args.type == 'subreads':
        row['mean_cov'] = int(gff_obj.getAttr(gff_id, "average_zmw_coverage"))
        row['retention_score'] = round(float(gff_obj.getValue(gff_id, 'score')), 4)
    lens = gff_obj.getAttr(gff_id, "IES_length")
    lens = [int(i) for i in lens.split("_")]
    if len(lens) == 1:
        row['length'] = lens[0]
    else:
        row['length'] = round(sum(lens)/ len(lens))

    if gff_obj.getAttr(gff_id, "ta_pointer_seq"):
        row['pointer'] = "ta"
    else:
        if gff_obj.getAttr(gff_id, "pointer_seq"):
            row['pointer'] = "pointer"
        else:
            row['pointer'] = "none"
    out.append(row)

# Dump JSON file with parsed data
with open(f"{args.out}.json", "w") as fh:
    fh.write(json.dumps(out, indent=4))

# Convert to Pandas data frame for easier plotting
df = pd.DataFrame.from_dict(out, dtype=object)

# IES retention score histogram
plt.figure(figsize=(8,6))
plt.hist([df.query('type == "ins"')['retention_score'],
          df.query('type == "del"')['retention_score']],
         label = ['ins','del'], stacked=True, bins=20)
plt.legend()
plt.xlabel("IES retention score")
plt.ylabel("Number of putative IESs")
plt.title("Retention score")
plt.savefig(f"{args.out}.ies_retention_score.{args.out_fmt}")

# IES lengths and pointer types
if args.hist_style == 'barstacked' or args.hist_style == 'bar':
    plt.figure(figsize=(8,6))
    plt.hist([df.query("pointer == 'ta'")["length"], 
              df.query("pointer == 'pointer'")["length"],
              df.query("pointer == 'none'")["length"]],
        label=['TA','other pointer','no pointer'], histtype=args.hist_style, bins=100)
    plt.legend()
    plt.title("IES length distribution")
    plt.savefig(f"{args.out}.ies_length_distribution.{args.out_fmt}")
# otherwise plot separately in facets
else:
    plt.figure(figsize=(8,18))
    plt.subplot(311)
    plt.hist(df.query("pointer == 'ta'")["length"],
             bins=100)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution - TA junction")

    plt.subplot(312)
    plt.hist(df.query("pointer == 'pointer'")["length"],
             bins=100)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution - other pointer")

    plt.subplot(313)
    plt.hist(df.query("pointer == 'none'")["length"],
             bins=100)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution - no pointer")

    plt.savefig(f"{args.out}.ies_length_distribution.{args.out_fmt}")

# IES lengths and pointer types - closeup
num_bins = args.hist_len_max - args.hist_len_min + 1

if args.hist_style == 'barstacked' or args.hist_style == 'bar':
    plt.figure(figsize=(8,6))
    plt.hist([df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'ta'")["length"], 
              df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'pointer'")["length"],
              df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'none'")["length"]],
        label=['TA','other pointer','no pointer'], histtype=args.hist_style, bins=num_bins)
    plt.legend()
    plt.title("IES length distribution (detail)")
    plt.savefig(f"{args.out}.ies_length_distribution_detail.{args.out_fmt}")
        
else:
    plt.figure(figsize=(8,18))
    plt.subplot(311)
    plt.hist(df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'ta'")["length"],
             bins=num_bins)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution (detail) - TA junction")

    plt.subplot(312)
    plt.hist(df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'pointer'")["length"],
             bins=num_bins)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution (detail) - other pointer")

    plt.subplot(313)
    plt.hist(df.query(f"length >= {args.hist_len_min} & length <= {args.hist_len_max} & pointer == 'none'")["length"],
             bins=num_bins)
    plt.xlabel("Length (bp)")
    plt.ylabel("Number of putative IESs")
    plt.title("IES length distribution (detail) - no pointer")

    plt.savefig(f"{args.out}.ies_length_distribution_detail.{args.out_fmt}")
