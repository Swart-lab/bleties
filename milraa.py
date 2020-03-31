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

import argparse
import sys
import pysam
from bleties import *

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
aln_format = "sam"
aln_mode = "r"
if args.sam:
    if args.bam:
        sys.exit("Error: Specify either SAM or BAM input, not both")
    else:
        aln_filename = args.sam
        aln_format = "sam"
        aln_mode = "r"
if args.bam:
    if args.sam:
        sys.exit ("Error: Specify either SAM or BAM input, not both")
    else:
        aln_filename=args.bam
        aln_format = "bam"
        aln_mode = "rb"


# Open SAM or BAM file 
alnfile = pysam.AlignmentFile(aln_filename, aln_mode)
# Initialize new IesRecords object to store putative IESs
iesrecords = IesRecords(alnfile, aln_format)

iesrecords.findPutativeIes(args.min_ies_length)

# print(iesrecords) # Print summary of IesRecords object
# iesrecords.dump() # Dump data to check

iesrecords.reportPutativeIes(args.min_break_coverage, args.out)

alnfile.close()
