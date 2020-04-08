#!/usr/bin/env python3

import argparse
import sys
from bleties import main

# Argument parser
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers()

# MILRAA
milraa_parser = subparsers.add_parser(name="milraa",
                                      description="MILRAA - Method of Identification by Long Read Alignment Anomalies")
milraa_parser.add_argument("--sam",
                           help="SAM file containing mapping, requires header")
milraa_parser.add_argument("--bam",
                           help="BAM file containing mapping, must be sorted and indexed")
milraa_parser.add_argument("--ref",
                           help="FASTA file containing genomic contigs used as reference for the mapping")
milraa_parser.add_argument("--out",
                           "-o",
                           nargs='?',
                           type=argparse.FileType("w"),
                           default=sys.stdout,
                           help="Path to write GFF3 file")
milraa_parser.add_argument("--out_fasta",
                            help="Path to write Fasta file of putative IES sequences")
milraa_parser.add_argument("--min_ies_length", # This parameter is hard-coded in the original ParTIES MIRAA
                            type=int,
                            default=25,
                            help="Minimum length of candidate IES")
milraa_parser.add_argument("--min_break_coverage", # For insertions
                            type=int,
                            default=10,
                            help="Minimum number of partially aligned reads to define a putative IES insertion breakpoint")
milraa_parser.add_argument("--min_del_coverage", # For deletions (sensu MILORD)
                            type=int,
                            default=10,
                            help="Minimum number of partially aligned reads to define a deletion relative to reference")
# milraa_parser.add_argument("--max_mismatch", # TODO: Not yet implemented
#                     type=int,
#                     default=10,
#                     help="Maximum mismatch in the alignment for a read to be used")
milraa_parser.add_argument("--dump",
                            action="store_true",
                            help="Dump contents of dict for troubleshooting")
# Assign function to this subparser
milraa_parser.set_defaults(func=main.milraa) 

# Parse arguments
args = parser.parse_args()
# Execute respective functions for each subparser
args.func(args)

"""The MIRAA module in ParTIES uses an alignment of Illumina reads vs somatic 
genome to look for breakpoints in read alignment. This script reimplements the 
MIRAA workflow for PacBio or other long read alignments.

Differences to Illumina alignments:
* One read may have multiple inserts
* Reads are not paired, insert size is not an issue
* Error rate of reads is expected to be higher
"""
