Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)
======================================================================

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for 
long read alignments. 

The required inputs are a ciliate MAC genome assembly, and a PacBio HiFi (high-
fidelity CCS reads) read library mapped onto that assembly. The mapper should
report a valid CIGAR string and NM tag (for number of mismatches) per aligned
read.


Usage
-----

Dependencies are specified as a Conda environment YAML file `env.yaml`.

For help, use the `-h` or `--help` option, with or without the subworkflow 
names:

```
./bleties.py --help
./bleties.py milraa --help
./bleties.py miser --help
./bleties.py milret --help
./bleties.py milcor --help
```

Scripts for plotting and visualizing data are in the `scripts/` subfolder.


Outline of workflow
-------------------

 * With MILRAA, identify putative IESs and IES junctions from mapping of PacBio
   HiFi/CCS reads to reference genome.
 * Screen for potential erroneous IES calls with MISER, curate the list of
   putative IES junctions.
 * Use curated IES junction list with MILRET to calculate IES retention scores.
 * Use curated IES junction list with MILCOR to calculate per-read IES retention
   scores and bin reads to MIC/MAC origin.


MILRAA - Method of Identification by Long Read Alignment Anomalies
------------------------------------------------------------------

Inputs:
 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed)
 * Reference assembly

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

Differences of long read to Illumina alignments:
 * One read may have multiple inserts -> Have to iterate through each alignment,
   rather than consider each only once
 * Reads are not paired, insert size is not an issue -> Ignore insert size
   parameter, expect reads to span entire IESs
 * Error rate of reads is expected to be higher -> Set higher `-max_mismatch`
   threshold

The output GFF3 file from MILRAA can be used to plot graphical summary of
predicted IESs with the script `scripts/milraa_gff_ies_plot.py`.


MISER - Method of IES Spurious or Erroneous Reporting
-----------------------------------------------------

Inputs:
 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed)
 * Reference assembly
 * Feature table of putative IES coordinates, produced by MILRAA

MISER takes an existing set of IES predictions, produced by MILRAA, and screens
it for potential mispredictions caused by paralogy, misassembly, or erroneous
mappings. 

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


MILRET - Method of IES Long-read Retention
------------------------------------------

Inputs:
 * Mapping of PacBio HiFi/CCS reads to reference genome
 * Reference assembly
 * Feature table (GFF3) of IES junctions, either from MILRAA output or third
   party tool

Reimplementation of MIRET in the ParTIES pipeline. MIRET takes the following
inputs: two reference genomes (somatic and germline), separate mapping files of
the same reads to the somatic and germline genomes, and coordinates of known
IES junctions. Reads that map with match to the somatic version are counted as
IES-, whereas reads that map with match to the germline version are counted as
IES+.

With long reads, we assume that IESs are spanned completely by most reads, vs
short reads, where reads are unlikely to completely span an IES insert. So we
do not count soft/hard clips on the ends of reads, only matches and inserts. We
also assume that the read mapper will handle mapping of reads containing inserts
properly. Therefore, we only map the reads to the somatic genome. Reads that do
contain IES+ forms will be reported as mappings with insert operations ("I" in
the CIGAR string). We then simply compare reads mapping with match to the
somatic genome at the IES junction (counted as IES-) to reads mapping with an
insert at the exact location of the IES junction (counted as IES+).

Differences of MILRET to MIRET:
 * MILRET uses only the somatic (IES-) as mapping reference, and assumes that 
   the mapper will correctly handle IES+ reads by mapping them with an insert 
   relative to the reference (I operation in the CIGAR string)
 * MILRET also counts and reports reads with additional deletions relative to
   the mapping reference. These may represent alternative excisions, 
   misassembly, or misalignments


MILCOR - Method of IES Long-read CORrelation
--------------------------------------------

Inputs:
 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed)
 * Feature table (GFF3) of IES junctions, either from MILRAA output or third
   party tool

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
