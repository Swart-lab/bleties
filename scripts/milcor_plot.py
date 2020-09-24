#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib as mpl
mpl.use("Agg") # allow run without X-server
import matplotlib.pyplot as plt


# Argument parser
parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--table",
        help="MILCOR output table")
parser.add_argument("--out", "-o", default="test",
        help="Prefix for output files")
args = parser.parse_args()


# Read table to Pandas object
df = pd.read_csv(args.table, sep="\t", header=0)
# header: ['qname', 'rname', 'start', 'end', 'ies_present', 'ies_absent']
df['qlen'] = df['end'] - df['start'] + 1 # coordinates are 1-based inclusive
df['retention_frac'] = df['ies_present'] / (df['ies_present'] + df['ies_absent'])


# Produce plot
plt.figure(figsize=(8,12))

plt.subplot(211)
plt.scatter(x=df['ies_present'],
        y=df['ies_absent'],
        alpha=0.2) # TODO: color by length of read
plt.xlabel("Number of IESs retained in read")
plt.ylabel("Number of IESs excised from read")
plt.title("IES retention per read")

plt.subplot(212)
plt.hist(df['retention_frac'])
plt.xlabel("Fraction of IESs retained per read")
plt.ylabel("Number of reads")
plt.title("IES retention per read as fraction of total per read")

plt.savefig(f"{args.out}.ies_retention_per_read.png")

