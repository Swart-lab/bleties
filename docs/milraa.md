MILRAA - Method of Identification by Long Read Alignment Anomalies
==================================================================

Reimplementation of MIRAA in the ParTIES pipeline. The original MIRAA is used
to identify IES retention in Paramecium genome by mapping Illumina reads
against a reference somatic (MAC) genome. It achieves this by an initial
Bowtie2 mapping step (local mode), then looks for breakpoints in the alignment,
which are identified as insert (I), soft clip (S), or hard clip (H) operations
in the CIGAR string of the alignment.

MILRAA reimplements the MIRAA concept for PacBio and other long read alignments.
The long reads must first be aligned to the reference genome (assumed to be
somatic possibly with some IES retentions that have assembled), assumption is
that the alignment is accurate.

The recommended aligner is [minimap2](https://github.com/lh3/minimap2), with
the options `--secondary=no --MD`, and with PacBio HiFi or CCS reads (use
option `-ax asm20`).


Input
-----
 * Ciliate MAC genome assembly, Fasta format.
 * PacBio HiFi or CCS read library mapping onto that assembly, sorted/indexed
   BAM format; mapper should report valid CIGAR string and NM tag.


Differences of long read to short read alignments
-------------------------------------------------

 * One read may have multiple inserts -> Have to iterate through each alignment,
   rather than consider each only once
 * Reads are not paired, insert size is not an issue -> Ignore insert size
   parameter, expect reads to span entire IESs
 * Error rate of reads is expected to be higher -> Set higher `-max_mismatch`
   threshold


Pointers and TA junctions
-------------------------

TK


Output
------

The output from MILRAA is a GFF3 file, with `MILRAA` in the `source` column
(column 2).

Putative IES junctions are reported where a minimum number of reads (threshold
defined by the `--min_break_coverage` option) contain an insert relative to the
reference assembly. In the GFF3 file, these are reported with
`internal_eliminated_sequence_junction` in the `type` column (column 3).
Junctions are features with zero length on the sequence; in GFF3 convention the
`end` coordinate (column 5) will equal the `start` coordinate (column 4), and
the junction is located to the right of that location coordinate.

Putative IES sequences are reported where a minimum number of reads (threshold
defined by the `--min_del_coverage` option) contain a deletion relative to the
reference assembly. In the GFF3 file, these are reported with
`internal_eliminated_sequence` in the `type` column.

The `score` column (column 6) reports the provisional IES retention score for
that putative IES feature. The IES retention score is defined as `IES+ / (IES+ +
IES-)` where `IES+` is the number of IES+ reads and `IES-` the number of IES-
reads for a given IES feature, in the mapped reads.

In addition, the `attributes` column (column 9) reports more information that is
relevant to IES analysis. Attributes are formatted as `tag=value` pairs,
delimited by semicolons `;`.

 * `ID` - Unique identifier for each IES reported
 * `IES_length` - Length of a given IES. If the `--fuzzy_ies` mode is used, and
     inserts of different lengths are used to define this IES, then the modal
     length is reported. If there is a tie between multiple values, then all
     tied values are reported, separated by underscore `_` character.
 * `cigar` - Summary of CIGAR operations that support a given reported IES. For
     example, if an IES junction is reported on the basis of 20 reads containing
     a 28 bp insert at that location, this is reported as `28I*20`. If
     `--fuzzy_ies` mode is used, and different operation lengths support a given
     reported IES, these are separated by space characters.
 * `average_coverage` - Total number of reads mapping to these coordinates. Only
     reported if input is a sorted, indexed BAM file.
 * `pointer_seq` - Possible pointer sequence, based on the mapped coordinates,
     if present. See above for discussion on pointer calling.
 * `ta_pointer_seq` - Pointer sequence starting with `TA`, if present.
 * `ta_pointer_start`, `ta_pointer_end` - Adjusted start and end coordinates
     relative to a TA boundary, if present.
 * `pp_pointer_seq` - Maximized pointer length, if present.
 * `pp_pointer_start`, `pp_pointer_end` - Adjusted start and end coordinates
     relative to a maximized pointer, if present.

The output GFF3 file from MILRAA can be used to plot graphical summary of
predicted IESs with the script `scripts/milraa_gff_ies_plot.py`.
