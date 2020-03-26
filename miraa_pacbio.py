#!/usr/bin/env python3

# MILRAA - Method of Identification by Long Read Alignment Anomalies
# The MIRAA module in ParTIES uses an alignment of Illumina reads
# vs somatic genome to look for breakpoints in read alignment
# This script reimplements the MIRAA workflow for PacBio or other
# long read alignments.
# Differences to Illumina alignments:
#  * One read may have multiple inserts
#  * Reads are not paired, insert size is not an issue
#  * Error rate of reads is expected to be higher

import re
import argparse
import json # for dumping data to debug
import sys
from collections import defaultdict

parser = argparse.ArgumentParser(description="MILRAA - MIRAA equivalent for long reads mappings, e.g. PacBio")
parser.add_argument("--sam",
                    nargs='?', 
                    type=argparse.FileType("r"), 
                    default=sys.stdin, 
                    help="SAM file containing mapping")
parser.add_argument("--out",
                    nargs='?',
                    type=argparse.FileType("w"),
                    default=sys.stdout,
                    help="Path to write GFF3 file")
parser.add_argument("--min_ies_length", # This parameter is hard-coded in the original ParTIES MIRAA
                    type=int,
                    default=25,
                    help="Minimum length of candidate IES")
parser.add_argument("--min_break_coverage",
                    type=int,
                    default=10,
                    help="Minimum number of partially aligned reads to define a breakpoint")
parser.add_argument("--max_mismatch", # TODO: Not yet implemented
                    type=int,
                    default=10,
                    help="Maximum mismatch in the alignment for a read to be used")
args = parser.parse_args()

# dict to store counts of detected inserts
# keys: contig -> startpos -> ies_length -> count
insDict = defaultdict(lambda: defaultdict( lambda: defaultdict (int)))

# Open SAM file and parse CIGAR string
for line in args.sam:
    linearr = line.split("\t") # Split line by tabs
    pos = int(linearr[3])
    rname = linearr[2]
    total_mismatch = 0
    nmmatch = re.search(r"NM:i:(\d+)", line)
    if nmmatch:
        total_mismatch = int(nmmatch.group(1))
    cigs = re.findall(r"\d+[\w\=]", linearr[5]) # Get CIGAR string
    # Initialize counters for how much query or reference seq is consumed
    # TODO: What about hard clipping?
    ref_consumed = 0
    que_consumed = 0
    total_i = 0 
    for cig in cigs: # for each cig in the string
        cigmatch = re.match(r"(\d+)([\w\=])",cig) # Get number and operation
        # We want to look for insert operations above a min length, 
        # map their positions on the reference, and count how many reads 
        # support a given putative insert
        # This requires that we count the operations that consume reference
        # and add it to the POS field
        if cigmatch.group(2) == "I": # If insert operation,
            ins_len = int(cigmatch.group(1)) # Length of the current insert
            total_i += ins_len # Add up the total insert length 
            if ins_len >= args.min_ies_length: # Check that insert is above min length
                # Get the start and end positions of the insert
                ins_pos_start = pos + ref_consumed # 1-based, insert is to the right of position
                insDict[rname][ins_pos_start][ins_len] += 2 # write to dict; counts as 2 because each insert has two ends
        # Count ref and query consumed _after_ the insert has been accounted for
        if cigmatch.group(2) in ['M', 'D', 'N', '=', 'X']:
            ref_consumed += int(cigmatch.group(1))
        if cigmatch.group(2) in ['M', 'I', 'S', '=', 'X']:
            que_consumed += int(cigmatch.group(1))
    if int(total_mismatch) - int(total_i) < 0: 
        # Sanity check - mismatches include inserts, but cannot be fewer than inserts
        print ("Uh-oh!")

# print(json.dumps(insDict, sort_keys=True, indent = 2))

# Parse the dict and report putative IESs above min coverage
for ctg in sorted(insDict):
    for ins_start in sorted(insDict[ctg]):
        for ins_len in sorted(insDict[ctg][ins_start]):
            if insDict[ctg][ins_start][ins_len] >= args.min_break_coverage:
                attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start),
                        "IES_length="+str(ins_len)
                       ] # TODO: cigar, average_coverage
                outarr = [str(ctg),        # 1 seqid
                          "MILRAA",        # 2 source
                          "segment",       # 3 type
                          str(ins_start),  # 4 start
                          str(ins_start),  # 5 end
                          str(insDict[ctg][ins_start][ins_len]), # 6 score - in this case, breakpoint counts
                          ".",             # 7 strand
                          ".",             # 8 phase
                          ";".join(attr)   # 9 attributes
                          ]
                args.out.write("\t".join(outarr))
