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
from bleties import getClips, getIndels

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
# keys: contig -> startpos -> endpos -> insert length -> evidence type -> count
insDict = defaultdict(
        lambda: defaultdict( 
            lambda: defaultdict ( 
                lambda: defaultdict(
                    lambda: defaultdict(int)
                    )
                )
            )
        )

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
    # total_i = 0 
    # Look for inserts that are on left or right ends of read (i.e. H or S clipping operations)
    cliparr = getClips(line.cigarstring, pos)
    if len(cliparr) > 0:
        for (clipstart, clipend, cliptype) in cliparr:
            insDict[rname][clipstart][clipend][0][cliptype] += 1
    # Look for inserts that are completely spanned by the read (i.e. I operations)
    indelarr = getIndels(line.cigarstring, pos, args.min_ies_length)
    if len(indelarr) > 0:
        for (indelstart, indelend, indellen, indeltype) in indelarr:
            insDict[rname][indelstart][indelend][indellen][indeltype] += 1
    # if int(total_mismatch) - int(total_i) < 0: 
        # Sanity check - mismatches include inserts, but cannot be fewer than inserts
        # print ("Uh-oh!")

# print(json.dumps(insDict, sort_keys=True, indent = 2))

# Parse the dict and report putative IESs above min coverage
# We only check breakpoints which are completely spanned by a read ("I" type)
# however we also report supporting counts from HSM and MSH type mappings
for ctg in sorted(insDict):
    for ins_start in sorted(insDict[ctg]):
        for ins_end in sorted(insDict[ctg][ins_start]):
            for ins_len in sorted(insDict[ctg][ins_start][ins_end]):
                for evidencetype in insDict[ctg][ins_start][ins_end][ins_len]:
                    countvalue = insDict[ctg][ins_start][ins_end][ins_len][evidencetype]
                    if evidencetype == "I":
                        if countvalue >= args.min_break_coverage:
                            # Prepare attributes list of key-value pairs
                            attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start)+"_"+str(ins_end),
                                    "IES_length="+str(ins_len)
                                   ]
                            # Extend the cigar evidence counts
                            attr.extend(["cigar="+cigartype+" "+str(insDict[ctg][ins_start][ins_end][ins_len][cigartype]) 
                                         for cigartype in sorted(insDict[ctg][ins_start][ins_end][ins_len])
                                        ])
                            # Get read coverage from BAM file; SAM does not allow random access
                            if aln_mode == "rb":
                                readcov = alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_end)) # TODO: Check for off-by-one errors
                                attr.append("average_coverage="+str(readcov))
                            outarr = [str(ctg),        # 1 seqid
                                      "MILRAA",        # 2 source
                                      "segment",       # 3 type
                                      str(ins_start),  # 4 start
                                      str(ins_end),  # 5 end
                                      str(countvalue), # 6 score - in this case, breakpoint counts for insert operation only
                                      ".",             # 7 strand
                                      ".",             # 8 phase
                                      ";".join(attr)   # 9 attributes
                                      ]
                            args.out.write("\t".join(outarr)+"\n")
                    elif evidencetype == "D":
                        if countvalue >= args.min_break_coverage:
                    # if insDict[ctg][ins_start][ins_end].get("D"): # If there is a deletion
                    #     if insDict[ctg][ins_start][ins_end]["D"] >= args.min_break_coverage:
                            del_len = int(ins_end) - int(ins_start) # TODO: Check for off-by-one errors
                            attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start)+"_"+str(ins_end),
                                    "IES_length="+str(del_len)
                                   ]
                            attr.append("cigar=D "+str(countvalue))
                            # Look for left- and rightclips that fall on the deletion boundaries
                            if insDict[ctg][ins_start][ins_start][ins_len].get("MHS") and int(insDict[ctg][ins_start][ins_start][ins_len]["MHS"]) > 0:
                                attr.append("cigar=MHS "+str(insDict[ctg][ins_start][ins_start][ins_len]["MHS"]))
                            if insDict[ctg][ins_end][ins_end][ins_len].get("HSM") and int(insDict[ctg][ins_end][ins_end][ins_len]["HSM"]) > 0:
                                attr.append("cigar=HSM "+str(insDict[ctg][ins_end][ins_end][ins_len]["HSM"]))
                            # Add coverage value
                            if aln_mode == "rb":
                                readcov = alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_end))
                                attr.append("average_coverage="+str(readcov))
                            outarr = [str(ctg),        # 1 seqid
                                      "MILRAA",        # 2 source
                                      "segment",       # 3 type
                                      str(ins_start),  # 4 start
                                      str(ins_end),  # 5 end
                                      str(countvalue), # 6 score - in this case, breakpoint counts for delete operation only
                                      ".",             # 7 strand
                                      ".",             # 8 phase
                                      ";".join(attr)   # 9 attributes
                                      ]
                            args.out.write("\t".join(outarr)+"\n")

alnfile.close()
