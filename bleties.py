#!/usr/bin/env python3

import argparse
import sys
import logging
from bleties import main

logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)
logging.info("Started BleTIES")
logging.info("Command line:")
logging.info(" ".join(sys.argv))
# Argument parser
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--log",
                    help="Log file",
                    default=sys.stderr)
subparsers = parser.add_subparsers()

# MILRAA -----------------------------------------------------------------------
"""The MIRAA module in ParTIES uses an alignment of Illumina reads vs somatic 
genome to look for breakpoints in read alignment. This script reimplements the 
MIRAA workflow for PacBio or other long read alignments.

Differences to Illumina alignments:
* One read may have multiple inserts
* Reads are not paired, insert size is not an issue
* Error rate of reads is expected to be higher
"""
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

# MILRET -----------------------------------------------------------------------
"""The MIRET pipeline in ParTIES compares mappings of the same reads to the 
somatic and germline genomes, at known IES junctions. Reads that map with match
to the somatic version are counted as IES-, reads that map with match to the 
germline version are counted as IES+. 

With long reads, we assume that IESs are spanned completely by most reads, vs
short reads, where reads are unlikely to completely span an IES insert. So we 
do not count soft/hard clips on the ends of reads, only matches and inserts. We
also assume that the read mapper will handle mapping of reads containing inserts
properly. Therefore, we only map the reads to the somatic genome. Reads that do
contain IES+ forms will be reported as mappings with insert operations ("I" in 
the CIGAR string). We then simply compare reads mapping with match to the 
somatic genome at the IES junction (counted as IES-) to reads mapping with an 
insert at the exact location of the IES junction (counted as IES+).

Points to note and address in the future:
 * Long reads are noisier and have many small indels. How do we distinguish 
   sequencing error from true inserts?
 * For same reason as above: What is minimum match length/quality before we 
   count a match?
 * Some junctions also exhibit deletions, which may be alternative excisions,
   misassembly, or misalignments. 
"""

milret_parser = subparsers.add_parser(name="milret",
                                      description="""MILRET - Method of IES Long-read RETention""")
milret_parser.add_argument("--bam",
                        help="BAM file containing mapping, must be sorted and indexed")
milret_parser.add_argument("--ref",
                        help="FASTA file containing genomic contigs used as reference for the mapping")
milret_parser.add_argument("--ies",
                        help="GFF3 file containing coordinates of IES junctions in MAC genome")
milret_parser.add_argument("--out",
                           "-o",
                           nargs='?',
                           type=argparse.FileType("w"),
                           default=sys.stdout,
                           help="Path to write table of retention scores per IES")
# Assign function to this subparser
milret_parser.set_defaults(func=main.milret)

# Parse arguments --------------------------------------------------------------
args = parser.parse_args()
# Execute respective functions for each subparser
args.func(args)
