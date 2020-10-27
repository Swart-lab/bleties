MILRET - Method of IES Long-read Retention
==========================================

The MILRET module calculates per-IES retention scores, either directly using the
predicted IES junctions predicted by [MILRAA](milraa.md), or IES junction
coordinates from third party software.


Inputs
------

 * Mapping of PacBio HiFi/CCS reads to reference genome
 * Reference assembly
 * Feature table (GFF3) of IES junctions, either from [MILRAA](milraa.md)
   output or third party tool


Differences of long read to short read alignments
-------------------------------------------------

MILRET is a reimplementation of the MIRET module in the ParTIES pipeline, that
was designed to work with short read data. MIRET takes the following inputs:
two reference genomes (somatic and germline), separate mapping files of the
same reads to the somatic and germline genomes, and coordinates of known IES
junctions. Reads that map with match to the somatic version are counted as
IES-, whereas reads that map with match to the germline version are counted as
IES+.

With long reads, we assume that IESs are spanned completely by most reads, vs
short reads, where reads are unlikely to completely span an IES insert. So we
do not count soft/hard clips on the ends of reads, only matches and inserts. We
also assume that the read mapper will handle mapping of reads containing
inserts properly, because most reads will contain predominantly MDS sequence
that will be correctly anchored in the somatic reference sequence. Therefore,
we only map the reads to the somatic genome. Reads that do contain IES+ forms
will be reported as mappings with insert operations (`I` in the CIGAR string).
We then simply compare reads mapping with match to the somatic genome at the
IES junction (counted as IES-) to reads mapping with an insert at the exact
location of the IES junction (counted as IES+).

 * MILRET uses only the somatic (IES-) as mapping reference, and assumes that 
   the mapper will correctly handle IES+ reads by mapping them with an insert 
   relative to the reference (I operation in the CIGAR string)
 * MILRET also counts and reports reads with additional deletions relative to
   the mapping reference. These may represent alternative excisions, 
   misassembly, or misalignments


Expected IES lengths
--------------------

Sequencing errors frequently cause short indels in reads, relative to the
reference assembly. Such erroneous indels may happen to overlap with IES
junction coordinates and be miscounted as evidence for an IES presence/absence.
When the `--use_ies_lengths` option is used, only indels with lengths matching
the expected IES length, as supplied in the input GFF3 file, will be counted as
a "true" IES retention. The expected IES length should be supplied with the
`IES_length` key in the `attributes` field of the GFF3 file, as done in the
output from [MILRAA](milraa.md).

Because of sequencing error, indel lengths may not match the expected IES length
exactly. This can be controlled with the `--length_threshold` option. For
example, a threshold of 0.05 means that sequences +/- 5% (inclusive) of the
expected length are still counted.


Output
------

The main output of MILRET is a table in TSV format `{OUT}.milret.tsv`, where
`{OUT}` is the output filename prefix supplied to the `--out` option. The TSV
table contains the following fields:

 * `ID` - Unique identifier for each IES junction, as supplied in the input GFF
   file.
 * `score` - Per-IES retention score, defined as IES+/(IES+ + IES-) read counts.
 * Remaining columns with names of CIGAR operations (e.g. `M`, `I`, `D`, `S`)
   report how many reads have that specific operation type at the coordinate
   corresponding to that IES junction. There may be a discrepancy between the
   reported counts and the retention score reported in the `score` column, if
   the option `--use_ies_lengths` is chosen. This is because a read may contain
   an indel operation at the correct coordinate, but the length does not match
   the expected IES length, e.g. in the case of indels caused by sequencing
   errors, which are typically much shorter than real IESs.

With the `--dump` option, internal data is dumped to JSON format for
troubleshooting to `{OUT}.milret.dump.json`
