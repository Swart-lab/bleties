MILTEL - Method of Long-read TELomere detection
===============================================

During development of MAC genomes in ciliates, chromosome breakage and addition
of telomeric repeats to the newly formed ends occurs. The extent of such
fragmentation varies between species, and in some species is known to be
regulated by conserved chromosome breakage sequences close to the breakage
sites.


Input data
----------

 * PacBio HiFi or CCS read library mapping onto MAC or MIC assembly,
   sorted/indexed BAM format; mapper should report valid CIGAR string and query
   sequences, and should perform soft clips instead of hard clips.
 * Known telomeric repeat unit sequence of the genome of interest, preferably in
   5' orientation.


Identification of telomere sequences from soft-clipped sequences
----------------------------------------------------------------

When a telomere-bearing sequence of MAC origin is mapped onto a MIC reference
sequence, the telomeric part does not align and is usually soft-clipped. Soft
clipping operations are limited to the ends of a read, i.e. there can be no more
than two soft clips on a single read.

MILTEL considers each mapped read with soft clips, and extracts the clipped
segment of the query sequence, the coordinates of the clip with respect to the
reference, and whether the unmapped sequence is to the left (5') or right (3')
of the reference coordinate.

Each clipped segment extracted above is searched for the user-provided telomeric
repeat sequence using
[NCRF](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder/), which can
find tandem repeats in the presence of noise from sequencing error. Where a
telomeric repeat (above a minimum length) is found, the gap distance from the
beginning/end of the telomeric repeat to the clipping junction is also counted
(in bp), as well as whether the telomere sequence is reverse complemented.


Output
------

MILTEL produces the following output files, where `{OUT}` is the output prefix
supplied to the `--out` option:

 * `{OUT}.miltel.telomeric.gff3` - Alternative telomere addition sites /
   chromosome breakage sites
 * `{OUT}.miltel.othercalls.gff3` - Other soft- or hard-clips that do not
   contain telomeric sequences
 * `{OUT}.miltel.othercalls.fasta` - Sequences of the above other clips

The GFF3 file for the alternative telomere addition sites has `MILTEL` in the
`source` column (column 2).  Coordinates where clipped sequence segments
contain telomeres are called putative `chromosome_breakage_site` in the `type`
column (column 3). Because this is a feature of zero length, the `start` and
`end` fields (columns 4 and 5) are equal, and the junction is to the right of
the coordinate, following GFF convention.

The `score` (column 6) reports the breakage score, which is the number of
telomere-bearing reads clipped at that specific coordinate, divided by the
total read coverage at that coordinate (as reported in the `average_coverage`
attribute, described below).

The `attributes` (column 9) contain the following fields:

 * `ID` - Unique identifier for each chromosome breakage site called.
 * `orientation` - Whether the clipped segment is to the `left` or to the
     `right` of the junction.
 * `telomere_sense` - Consensus (majority rule) on whether the clipped telomere
     sequences are sense (`+`) or reverse complement (`-`) of the supplied
     repeat unit sequence. If the supplied repeat sequence is the 5' telomere
     sequence then the `left` orientation should correspond to the `+` telomere
     sense, and `right` to `-`.
 * `telomere_gap_average` - Mean gap distance between beginning of telomeric
     repeats and the clipping junction.
 * `telomere_senses` - String reporting the telomere sense orientations when
     there is more than one read clipped at that coordinate. For example, if
     there are 2 reads with `+` and 1 read with `-`, then the string is
     `telomere_senses=2*+ 1*-`.
 * `telomere_gaps` - String reporting the gap distances for each clipped read.
     For example if there were 2 reads with gap 0 and one with gap 4, then the
     string is `telomere_gaps=0*2 4*1`
 * `average_coverage` - Total read coverage at the coordinate, excluding
     secondary and supplementary mappings. Calculated the same way as the
     corresponding field in MILRAA output.

With the `--dump` option, internal data are dumped in JSON format for
troubleshooting to: `{OUT}.miltel.dump.json`.
