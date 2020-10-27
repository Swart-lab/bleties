MILCOR - Method of IES Long-read CORrelation
============================================

With long reads (>1 kbp) it is possible to count IES retention at the level of
individual reads. In PacBio or Nanopore sequencing, libraries are prepared
without amplification (e.g. by PCR), so the reads represent original molecules,
without the possibility of PCR chimerism. We could therefore potentially
classify reads into MIC-origin or MAC-origin, in the case of vegetative cells,
or examine the dynamics of IES excision in developing MACs.

MILCOR reports a per-read IES "retention" score that complements the per-IES
retention score reported by [MILRET](milret.md). This is not possible with
short read sequencing where reads typically do not span an entire IES. In the
calculation of the per-IES "retention" score, reads that do not span at least
one defined IES junction site are not counted.


Inputs
------

 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed).
 * Feature table (GFF3) of IES junctions, either from [MILRAA](milraa.md)
   output or third party tool.

<!--
Parameters
---------

-->

Terminology
-----------

The IES retention score was originally defined as a measure of IES excision
efficiency in a population of cells derived from an experimental knockdown of a
candidate gene involved in IES excision. 100% excision efficiency would lead to
retention score of 0.

In the context of a single read, the term "retention" is arguably
inappropriate, because it was originally defined for a population of sequences,
rather than for a single sequence which may originate from a MAC, MIC, or
developing MAC. Nonetheless, the per-read measure will still be referred to as
a "retention score" here for convenience and to draw parallels to the per-IES
retention score.


Read binning
------------

MAC-derived sequences are expected to have a per-read retention score of 0,
whereas MIC-derived sequences should have a per-read retention score of 1. It
is therefore possible to bin reads into putatively MAC or MIC origin, based on
retention scores of close to 0 or 1 respectively. This is done with the `--bin`
option. The `--bin_threshold` option sets the minimum excision/retention
required. For example, `--bin_threshold 0.9` (the default) means that sequences
with per-read score $\geq 0.9$ will be binned as MIC, and score $\leq 0.1$ will
be binned as MAC.


Output
------

The main output from MILCOR is a table in TSV format, with IES presence/absence
statistics per read, with filename `{OUT}.milcor.tsv` where `{OUT}` is the
output filename prefix supplied to the `--out` option. The table has the
following fields:

 * `qname` - Name of the read, from BAM file
 * `rname` - Name of contig/scaffold in reference with the primary mapping of
   this read.
 * `start`, `end` - Coordinates on the reference where the read maps, from BAM
   file
 * `ies_present` - Number of IES junctions within those reference coordinates
   where the read contains an insert (i.e. IES not excised)
 * `ies_absent` - Number of IES junctions within those reference coordinates
   where the read does not contain an insert (i.e. IES excised)

The TSV file can be used as input to plot a graphical summary of the per-read
IES retention scores, with the script `milcor_plot.py`.

With the `--dump` option, internal data is dumped to JSON format for
troubleshooting to `{OUT}.milcor.dump.json`.

With the `--bin` option, binned reads are reported in Fasta format to the
following files (see "Read binning" above):

 * `{OUT}.milcor_bin_MAC.fasta` - IES- reads likely to be of MAC origin, with
   per-read retention score below $1 - b$ where $b$ is the threshold value
   supplied to `--bin_threshold`
 * `{OUT}.milcor_bin_MIC.fasta` - IES+ reads likely of MIC origin
 * `{OUT}.milcor_bin_other.fasta` - Reads with per-read retention score below
   the thresholds for either MIC or MAC bins.
 * `{OUT}.milcor_bin_noies.fasta` - Reads that do not span any annotated IES
   junctions on their mapping to the reference genome, and hence cannot be
   placed into either MIC or MAC bins.
