#!/usr/bin/env python3

from bleties import SharedFunctions
from sys import argv
import re
import json

gff_file = argv[1]

gff_obj = SharedFunctions.Gff()
gff_obj.file2gff(gff_file)

out = []
for gff_id in gff_obj:
    row = {}
    row['id'] = gff_id
    if gff_obj.getValue(gff_id, "type") == "junction":
        row['type'] = "ins"
    elif gff_obj.getValue(gff_id, "type") == "region":
        row['type'] = "del"
    row['mean_cov'] = int(gff_obj.getAttr(gff_id, "average_coverage"))
    cigar = gff_obj.getAttr(gff_id, "cigar")
    cigar_covs = re.findall(r"\d+[ID]\*(\d+)", cigar)
    cigar_total = sum([int(i) for i in cigar_covs])
    if row['type'] == "ins":
        iesplus = cigar_total
    elif row['type'] == "del":
        iesplus = row['mean_cov'] - cigar_total
    row['retention_score'] = iesplus / row['mean_cov']
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

print(json.dumps(out, indent=4))
