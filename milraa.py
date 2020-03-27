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
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser(description="MILRAA - MIRAA equivalent for long reads mappings, e.g. PacBio")
parser.add_argument("--sam",
                    help="SAM file containing mapping, requires header")
parser.add_argument("--bam",
                    help="BAM file containing mapping, must be sorted and indexed")
parser.add_argument("--out",
                    "-o",
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

# Check that only either SAM or BAM specified
aln_filename = "-"
aln_mode = "r"
if args.sam:
    if args.bam:
        sys.exit("Error: Specify either SAM or BAM input, not both")
    else:
        aln_filename = args.sam
        aln_mode = "r"
if args.bam:
    if args.sam:
        sys.exit ("Error: Specify either SAM or BAM input, not both")
    else:
        aln_filename=args.bam
        aln_mode = "rb"

# dict to store counts of detected inserts keyed by evidence type
# keys: contig -> startpos -> evidence type -> count
insDict = defaultdict(lambda: defaultdict( lambda: defaultdict (int)))
# dcit to store counts of detected inserts keyed by insert length
# keys: contig -> startpos -> ies_length -> count
insLenDict = defaultdict(lambda: defaultdict( lambda: defaultdict (int)))

# Open SAM file and parse CIGAR string
alnfile = pysam.AlignmentFile(aln_filename, aln_mode)
for line in alnfile:
    # linearr = line.split("\t") # Split line by tabs
    pos = int(line.reference_start) + 1 # Convert from 0-based numbering in pysam to 1-based in GFF3 and SAM
    rname = line.reference_name # Get reference name
    total_mismatch = line.get_tag("NM") # Get number of mismatches
    cigs = re.findall(r"\d+[\w\=]", line.cigarstring) # Get CIGAR string
    # Initialize counters for how much query or reference seq is consumed
    ref_consumed = 0
    que_consumed = 0
    total_i = 0 
    # Look for inserts that are on left or right ends of read (i.e. H or S clipping operations)
    leftclipmatch = re.match(r"(\d+)[HS](\d+)M", line.cigarstring)
    rightclipmatch = re.search(r"(\d+)M(\d+)[HS]$", line.cigarstring)
    if leftclipmatch:
        insDict[rname][pos]['HSM'] += 1
    if rightclipmatch:
        refconsumematches = re.findall(r"(\d+)[MDN\=X]", line.cigarstring) # find all ref-consuming operations
        rightclip_refconsumed = sum([int(refconsumematch) for refconsumematch in refconsumematches])
        insDict[rname][pos + rightclip_refconsumed]['MHS'] += 1
    # Look for inserts that are completely spanned by the read (i.e. I operations)
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
                insLenDict[rname][ins_pos_start][ins_len] += 1 # record the putative IES length spanned by the read
                insDict[rname][ins_pos_start]["I"] += 2 # record as "I" type count, like MIRAA; counts as 2 because each insert has two ends
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
# We only check breakpoints which are completely spanned by a read ("I" type)
# however we also report supporting counts from HSM and MSH type mappings
for ctg in sorted(insLenDict):
    for ins_start in sorted(insLenDict[ctg]):
        for ins_len in sorted(insLenDict[ctg][ins_start]):
            if insLenDict[ctg][ins_start][ins_len] >= args.min_break_coverage:
                # Prepare attributes list of key-value pairs
                attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start),
                        "IES_length="+str(ins_len)
                       ]
                # Extend the cigar evidence counts
                attr.extend(["cigar="+cigartype+" "+str(insDict[ctg][ins_start][cigartype]) for cigartype in sorted(insDict[ctg][ins_start])])
                # Get read coverage from BAM file; SAM does not allow random access
                if aln_mode == "rb":
                    readcov = alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_start))
                    attr.append("average_coverage="+str(readcov))
                outarr = [str(ctg),        # 1 seqid
                          "MILRAA",        # 2 source
                          "segment",       # 3 type
                          str(ins_start),  # 4 start
                          str(ins_start),  # 5 end
                          str(insLenDict[ctg][ins_start][ins_len]), # 6 score - in this case, breakpoint counts
                          ".",             # 7 strand
                          ".",             # 8 phase
                          ";".join(attr)   # 9 attributes
                          ]
                args.out.write("\t".join(outarr)+"\n")

alnfile.close()
