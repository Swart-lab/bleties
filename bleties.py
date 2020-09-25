#!/usr/bin/env python3

import argparse
import sys
import logging

from bleties import main

# Argument parser
parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--log", default="bleties.log",
    help="Path to write log file")
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
    description="""
    MILRAA - Method of Identification by Long Read Alignment Anomalies
    """,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Input arguments
milraa_parser.add_argument("--sam",
    help="SAM file containing mapping, requires header")
milraa_parser.add_argument("--bam",
    help="BAM file containing mapping, must be sorted and indexed")
milraa_parser.add_argument("--ref",
    help="""
    FASTA file containing genomic contigs used as reference for the mapping
    """)

# Output arguments
milraa_parser.add_argument("--out",
    "-o",
    default="milraa.test",
    help="Output filename prefix")

# Params for reporting candidate IESs
milraa_parser.add_argument("--junction_flank",
    type=int,
    default=5,
    help="Length of flanking sequence to report to --out_junction")
milraa_parser.add_argument("--min_ies_length", # This parameter is hard-coded in the original ParTIES MIRAA
    type=int,
    default=25,
    help="Minimum length of candidate IES")
milraa_parser.add_argument("--min_break_coverage", # For insertions
    type=int,
    default=10,
    help="""
    Minimum number of partially aligned reads to define a putative IES insertion
    breakpoint
    """)
milraa_parser.add_argument("--min_del_coverage", # For deletions (sensu MILORD)
    type=int,
    default=10,
    help="""
    Minimum number of partially aligned reads to define a deletion relative to
    reference
    """)
# milraa_parser.add_argument("--max_mismatch", # TODO: Not yet implemented
#    type=int,
#    default=10,
#    help="Maximum mismatch in the alignment for a read to be used")

# Params for fuzzy length reporting candidate IESs
milraa_parser.add_argument("--fuzzy_ies",
    action="store_true",
    help="""
    Allow lengths of inserts to differ slightly when defining putative IES,
    otherwise insert lengths must be exactly the same.
    """)
milraa_parser.add_argument("--cluster_dist",
    default=0.05,
    type=float,
    help="""
    Sequence identity distance limit for clustering putative IESs together. 
    Recommended settings: 0.05 for PacBio HiFi reads, ??? for PacBio CLR reads.
    Not yet tested extensively.
    """)

# Others
milraa_parser.add_argument("--dump",
    action="store_true",
    help="Dump contents of dict for troubleshooting")

# Assign function to this subparser
milraa_parser.set_defaults(func=main.milraa)

# MISER ------------------------------------------------------------------------
"""MISER takes an existing set of IES predictions, produced by MILRAA, and
screens it for potential mispredictions caused by paralogy, misassembly, or
erroneous mappings.

For each putative IES (insertion or deletion), the set of reads mapping to that
site is found, and split into two subsets: those containing the indel and those
without. For each subset, the mean percent mismatch of alignments vs. the
reference is taken.

 * If either subset has high (>5%) mismatch rate, "high error" is reported.
 * If the subset with indel has a significantly higher mismatch than the subset
   without, possible paralog is reported.
 * If the subset without indel has a significantly higher mismatch, possible
   misassembly is reported.
 * Otherwise the putative IES is "ok".
"""
miser_parser = subparsers.add_parser(name="miser",
    description="MISER - Method of IES Spurious or Erroneous Reporting",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Input arguments
miser_parser.add_argument("--sam",
    help="SAM file containing mapping, requires header")
miser_parser.add_argument("--bam",
    help="BAM file containing mapping, must be sorted and indexed")
miser_parser.add_argument("--ref",
    help="""
    FASTA file containing genomic contigs used as reference for the mapping""")
miser_parser.add_argument("--gff",
    help="GFF file containing coordinates for putative IESs")

# Output arguments
miser_parser.add_argument("--out",
    "-o",
    nargs='?',
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="""
    Path to write report statistics on possibly spurious IESs due to misassembly
    or mapped paralogs, defaults to STDOUT
    """)
miser_parser.add_argument("--split_gff",
    action="store_true",
    help="""
    Split input GFF entries into separate files for each category (ok,
    misassembly, paralog, ...), using input GFF filename as prefix
    """)

# Test parameters
miser_parser.add_argument("--spurious_ies_test",
    type=str,
    default="mann-whitney",
    help="""
    Test to use to evaluate spurious IESs by mismatch percentage comparisons,
    either \"mann-whitney\" (Mann-Whitney's U) or \"t\" (Ward's t-test)
    """)
miser_parser.add_argument("--spurious_ies_pvalue",
    type=float,
    default=0.05,
    help="""
    P-value cutoff (uncorrected) to use for spurious IES mismatch test; the
    Bonferroni correction will be applied depending on the number of tests
    (number of putative IESs) performed
    """)

# Assign function to this subparser
miser_parser.set_defaults(func=main.miser)

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
 * Deletions are not counted because the reference is assumed to be IES-free.
   Entries in the input GFF that define regions (i.e. not junctions) will be
   ignored.
"""

milret_parser = subparsers.add_parser(name="milret",
    description="""
    MILRET - Method of IES Long-read RETention
    """,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Input arguments
milret_parser.add_argument("--bam",
    help="BAM file containing mapping, must be sorted and indexed")
milret_parser.add_argument("--ref",
    help="""
    FASTA file containing genomic contigs used as reference for the mapping""")
milret_parser.add_argument("--ies",
    help="GFF3 file containing coordinates of IES junctions in MAC genome")
milret_parser.add_argument("--use_ies_lengths", action="store_true",
    help="""
    Only count inserts that match IES lengths reported in the input GFF file.
    This assumes that the input GFF file is produced by BleTIES MILRAA""")
milret_parser.add_argument("--length_threshold", type=float, default=0.05,
    help="""
    Length threshold to count matching IES length, if option --use_ies_lengths
    is applied""")

# Output arguments
milret_parser.add_argument("--out",
    "-o",
    default="milret.test",
    help="Path to write table of retention scores per IES")
milret_parser.add_argument("--dump", action="store_true",
    help="Dump contents of retention score objects to JSON file, for troubleshooting")

# Assign function to this subparser
milret_parser.set_defaults(func=main.milret)


# MILCOR -----------------------------------------------------------------------
"""MILCOR - Method of IES Long-read CORrelation
With long reads (>1 kbp) it is possible to count IES retention at the level of
individual reads. In PacBio or Nanopore sequencing, libraries are prepared
without amplification (e.g. by PCR), so the reads represent original molecules,
without the possibility of PCR chimerism. We could therefore potentially
classify reads into MIC-origin or MAC-origin, in the case of vegetative cells,
or examine the dynamics of IES excision in developing MACs.

MILCOR reports a per-read IES retention score that complements the per-IES
retention score reported by MILRET. This is not possible with short read
sequencing where reads typically do not span an entire IES. In the calculation
of the per-IES retention score, reads that do not span at least one defined IES
junction site are not counted.
"""

milcor_parser = subparsers.add_parser(name="milcor",
        description="MILCOR - Method of IES Long-read CORrelation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Input arguments
milcor_parser.add_argument("--bam",
    help="BAM file containing mapping, must be sorted and indexed")
milcor_parser.add_argument("--ies",
    help="GFF3 file containing coordinates of IES junctions in MAC genome")

# Output arguments
milcor_parser.add_argument("--out",
    "-o",
    default="milcor.test",
    help="Path to write table of per-read IES retention scores")
milcor_parser.add_argument("--dump", 
    action="store_true",
    help="Dump contents of IES correlation objects to JSON file, for troubleshooting")

milcor_parser.set_defaults(func=main.milcor)

# Parse arguments --------------------------------------------------------------
args = parser.parse_args()

# Logging
logging.basicConfig(format='[%(asctime)s] %(name)-12s %(levelname)-8s %(message)s',
        level=logging.DEBUG,
        filename=args.log,
        filemode="a")
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)-24s %(levelname)-8s %(message)s")
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# Execute respective functions for each subparser
args.func(args)
