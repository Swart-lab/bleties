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

# Relative import relative to `scripts` folder
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from bleties import SharedFunctions


# Argument parser
parser = argparse.ArgumentParser(
        description="Plot IES statistics from MILRAA output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gff",
        help="GFF file produced by BleTIES MILRAA module")
parser.add_argument("--out", "-o", default="test",
        help="Prefix for output files")
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
    row['mean_cov'] = int(gff_obj.getAttr(gff_id, "average_coverage"))
    cigar = gff_obj.getAttr(gff_id, "cigar")
    cigar_covs = re.findall(r"\d+[ID]\*(\d+)", cigar)
    cigar_total = sum([int(i) for i in cigar_covs])
    if row['type'] == "ins":
        iesplus = cigar_total
    elif row['type'] == "del":
        iesplus = row['mean_cov'] - cigar_total
    row['retention_score'] = round(float(iesplus / row['mean_cov']), 4)
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
plt.savefig(f"{args.out}.ies_retention_score.png")

# IES lengths and pointer types
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

plt.savefig(f"{args.out}.ies_length_distribution.png")


# IES lengths and pointer types - closeup
plt.figure(figsize=(8,18))
plt.subplot(311)
plt.hist(df.query("length > 25 & length <= 400 & pointer == 'ta'")["length"],
         bins=375)
plt.xlabel("Length (bp)")
plt.ylabel("Number of putative IESs")
plt.title("IES length distribution (detail) - TA junction")
plt.savefig(f"{args.out}.ies_stats_plots.png")

plt.subplot(312)
plt.hist(df.query("length > 25 & length <= 400 & pointer == 'pointer'")["length"],
         bins=375)
plt.xlabel("Length (bp)")
plt.ylabel("Number of putative IESs")
plt.title("IES length distribution (detail) - other pointer")
plt.savefig(f"{args.out}.ies_stats_plots.png")

plt.subplot(313)
plt.hist(df.query("length > 25 & length <= 400 & pointer == 'none'")["length"],
         bins=375)
plt.xlabel("Length (bp)")
plt.ylabel("Number of putative IESs")
plt.title("IES length distribution (detail) - no pointer")

plt.savefig(f"{args.out}.ies_length_distribution_detail.png")
