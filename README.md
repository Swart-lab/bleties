# Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)

This is a reimplementation of ParTIES for long read alignments. 

## Usage

Dependencies are specified as a Conda environment YAML file `env.yaml`.

For help, use the `-h` or `--help` option, with or without the subworkflow 
names:

```
./bleties.py --help
./bleties.py milraa --help
./bleties.py milret --help
```

## MILRAA - Method of Identification by Long Read Alignment Anomalies

Reimplementation of MIRAA in the ParTIES pipeline. MIRAA is used to identify IES
retention in Paramecium genome by mapping Illumina reads against a reference
somatic (MAC) genome. It achieves this by an initial Bowtie2 mapping step (local
mode), then looks for breakpoints in the alignment, which are identified as
insert (I), soft clip (S), or hard clip (H) operations in the CIGAR string of
the alignment.

MILRAA reimplements the MIRAA concept for PacBio and other long read alignments.
The long reads must first be aligned to the reference genome (assumed to be
somatic possibly with some IES retentions that have assembled), assumption is
that the alignment is accurate.

Differences of long read to Illumina alignments:
 * One read may have multiple inserts -> Have to iterate through each alignment,
   rather than consider each only once
 * Reads are not paired, insert size is not an issue -> Ignore insert size
   parameter, expect reads to span entire IESs
 * Error rate of reads is expected to be higher -> Set higher `-max_mismatch`
   threshold

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

