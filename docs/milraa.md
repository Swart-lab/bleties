MILRAA - Method of Identification by Long Read Alignment Anomalies
==================================================================

The MILRAA module predicts likely IES insertion/deletion coordinates from a
mapping of long sequencing reads onto a reference genome assembly. The reads
must first be aligned to the reference, which is assumed to be somatic (MAC) but
possibly with some IESs retained in the assembly. The alignment is assumed to be
accurate. The outputs are the coordinates of the predicted IES indel positions,
and IES consensus sequences.

The recommended aligner is [minimap2](https://github.com/lh3/minimap2), with
the options `--secondary=no --MD`, and with PacBio HiFi or CCS reads (use
option `-ax asm20`).


Input data
----------

 * Ciliate MAC genome assembly, Fasta format.
 * PacBio HiFi or CCS read library mapping onto that assembly, sorted/indexed
   BAM format; mapper should report valid CIGAR string and NM tag.

<!--
Arguments
---------

### Inputs

### Outputs

### Reporting parameters

### Other
-->

Differences of long read to short read alignments
-------------------------------------------------

MILRAA is a reimplementation of the MIRAA module of the original ParTIES
pipeline. The original MIRAA was used for Illumina short read data, whereas
MILRAA should be used with long read data. ParTIES MIRAA uses an initial
Bowtie2 mapping step (local mode), then looks for breakpoints in the alignment,
which are identified as insert (I), soft clip (S), or hard clip (H) operations
in the CIGAR string of the alignment.

The insert size for Illumina paired end read libraries is typically around 400
bp, so on average each read-pair would span at most one IES. Most read pairs are
unlikely to span a complete IES nor contain the entire IES sequence. In
comparison, a PacBio long read with lengths 10 kbp or more may span multiple
IESs, and the entire IES sequence can potentially be read out form the sequence.
Therefore, MILRAA has the following differences to MIRAA:

 * Each read potentially covers multiple indels, unlike MIRAA which considers
   one read as support for at most one IES.
 * Reads are not paired, so "insert size" is not an issue, and the entire insert
   sequence can be extracted from the reads.
 * However, the error rate of PacBio reads, even HiFi/CCS reads, is expected to
   be higher, so extracted insert sequences are aligned to get a consensus IES
   sequence, and a minimum indel length is required for IES calling (because
   many short indels are expected from sequencing error).


Pointers and TA junctions
-------------------------

IESs can be bounded by TA junctions or pointers (tandem repeats). One copy of
the pointer or TA is excised with the IES, while another copy remains on the
macronuclear destined sequence (MDS).

MILRAA uses the following conventions for coordinates, although they may not
correspond to the way that the physical bases are excised in the actual DNA
molecule:

 * The coordinate of an IES junction (i.e. IES has been excised from the
   reference genome) is to the left of a pointer or TA sequence in the
   reference genome.
 * The coordinate of an IES region (i.e. IES is retained in the reference
   genome) is to the left of the first pointer or TA sequence. The leftmost
   pointer/TA is considered the IES copy, and the rightmost pointer/TA as the
   MDS copy.

The coordinates reported in the `start` and `end` columns of the MILRAA GFF3
output are the coordinates based on the insert/deletion operations reported by
the mapper. However, if an indel is flanked by tandem repeats, there is
ambiguity as to where the "true" indel is, because it is not possible to tell
from the sequence alone where in the repeated sequence the physical cut
occurred.

For example, consider the following sequence containing an IES (upper case)
surrounded by MDS (lower case), the junctions are indicated with the pipe
character `|`:

```
gcgc|TAATGGTGCC|taatccgc (MDS + IES)

gcgc|----------|taatccgc (MDS only)

----|TAATGGTGCC|-------- (IES inferred sequence)
```

Notice the pointer sequence `TAAT` repeated in the left side in the IES, and to
the right of the junction in the MDS.

The first sequence (MDS + IES), when mapped onto the MDS sequence, would have the CIGAR
string `4M10I8M`. 

However, it is also possible to map it another way, with the CIGAR string
`8M10I4M`. The `TAAT` sequence now lies on the right side in the IES, and to
the left of the junction in the MDS. The MDS + IES and MDS-only sequences have
not changed, only our decision of where to place the junction of the inferred
indel.

```
gcgctaat|GGTGCCTAAT|ccgc (MDS + IES)

gcgctaat|----------|ccgc (MDS only)

--------|GGTGCCTAAT|---- (IES inferred sequence)
```

Because the insert is flanked by the repeated pointer sequence, alternative
mappings of `5M10I7M`, `6M10I6M`, and `7M10I5M` are also possible. To deal with
this ambiguity, MILRAA does the following:

 * If the putative IES indel, as reported by the mapper's original output, is
   flanked by repeats, these are reported as putative pointer sequences, along
   with the original indel coordinates from the mapper.
 * If the mapper reports the indel with the pointer/TA on the right side in the
   IES region, then the coordinates are adjusted so that it lies in the left
   side in the IES.
 * If the putative pointer contains a "TA" sequence, this is reported as a
   possible TA junction, and the adjusted coordinates for this TA junction are
   reported.
 * If it is possible to extend the repeats on their left sides to be longer
   than what the mapper reports, then the longest possible repeat (pointer) is
   reported, as well as the corresponding adjusted coordinates.

The pointer sequences and adjusted coordinates are reported in the `attributes`
field of the MILRAA GFF3 file, as described below.


Fuzzy IES length matching
-------------------------

The default mode is "strict", i.e. to define a given IES junction, all the
insert operations at a given coordinate must have the same length. However,
this does not account for sequencing errors that may introduce indels, nor the
possibility of different excised sequences at the same coordinate.

With the `--fuzzy_ies` option, insert sequences (putative IES sequences) at a
given coordinate are extracted from the mapped reads. These are pairwise
aligned and clustered by sequence identity. The clusters are then split at a
identity threshold given by the `--cluster_dist` option (default 0.05, i.e.
5%). Clusters are then used to define an IES junction if they meet the minimum
coverage threshold. The consensus IES sequence is produced from a multiple
sequence alignment (with Muscle) of the cluster members. Gaps are reported in
the output Fasta file (with the `-` character), but gaps are ignored when
predicting a pointer sequence.


Output
------

MILRAA produces the following output files by default (`{OUT}` is the output
filename prefix, supplied with the `--out` option):

 * `{OUT}.milraa_ies.gff3`
 * `{OUT}.milraa_ies.fasta`

With the `--fuzzy_ies` option, the following additional outputs are produced:

 * `{OUT}.milraa_ies_fuzzy.gff3`
 * `{OUT}.milraa_ies_fuzzy.fasta`


### MILRAA GFF3 output

The main output from MILRAA is a GFF3 file, with `MILRAA` in the `source` column
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


### MILRAA Fasta output

The Fasta file contains the consensus IES sequences, with headers corresponding
to the `ID` field in the GFF3 file.
