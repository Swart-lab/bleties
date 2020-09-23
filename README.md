# Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for 
long read alignments. 

The required inputs are a ciliate MAC genome assembly, and a PacBio HiFi (high-
fidelity CCS reads) read library mapping onto that assembly. The mapper should
report valid CIGAR string and NM tag (for number of mismatches) per aligned
read.

## Usage

Dependencies are specified as a Conda environment YAML file `env.yaml`.

For help, use the `-h` or `--help` option, with or without the subworkflow 
names:

```
./bleties.py --help
./bleties.py milraa --help
./bleties.py miser --help
./bleties.py milret --help
```

Utility scripts are in the `scripts/` subfolder.

## MILRAA - Method of Identification by Long Read Alignment Anomalies

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

## MISER - Method of IES Spurious or Erroneous Reporting

MISER takes an existing set of IES predictions, produced by MILRAA, and 
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

## MILRET - Method of IES Long-read Retention

Reimplementation of MIRET in the ParTIES pipeline. MIRET takes two reference
genomes: somatic and germline, and separate mapping files of the same reads to
the somatic and germline genomes, to identify reads that map to the IES+
(germline) vs. those that map to the IES- (somatic) forms. 

Differences of MILRET to MIRET:
 * MILRET uses only the somatic (IES-) as mapping reference, and assumes that 
   the mapper will correctly handle IES+ reads by mapping them with an insert 
   relative to the reference (I operation in the CIGAR string)
 * MILRET also counts and reports reads with additional deletions relative to
   the mapping reference. These may represent alternative excisions, 
   misassembly, or misalignments

