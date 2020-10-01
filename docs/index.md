Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)
======================================================================

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for 
long read alignments. 


Input data
----------

 * Ciliate MAC genome assembly, Fasta format.
 * PacBio HiFi or CCS read library mapping onto that assembly, sorted/indexed
   BAM format; mapper should report valid CIGAR string and NM tag.


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

Refer to the individual module pages for further information.

 * [MILRAA](milraa.md) -- Identify putative IESs and IES junctions from a
   mapping of PacBio HiFi/CCS reads to reference genome assembly.
 * [MISER](miser.md) -- Screen for potentially erroneous IES calls, to curate
   the list of putative IES junctions.
 * [MILRET](milret.md) -- Use curated IES junctions, or IES junction coordinates
   from third-party tools, to calculate IES retention scores from mapping.
 * [MILCOR](milcor.md) -- Calculate per-read IES retention scores and bin reads
   as MIC/MAC in origin.

