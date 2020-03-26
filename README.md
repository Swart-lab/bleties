# Reimplementation of MIRAA for long read alignments

MIRAA is a module in the ParTIES pipeline, used to identify IES retention in Paramecium genome by mapping Illumina reads against a reference somatic (MAC) genome. It achieves this by an initial Bowtie2 mapping step (local mode), then looks for breakpoints in the alignment, which are identified as insert (I), soft clip (S), or hard clip (H) operations in the CIGAR string of the alignment.

This script reimplements the MIRAA concept for PacBio and other long read alignments. 

Differences of long read to Illumina alignments:
 * One read may have multiple inserts
 * Reads are not paired, insert size is not an issue
 * Error rate of reads is expected to be higher
