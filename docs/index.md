Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)
======================================================================

![BLETIES logo](./bleties_logo.png)

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for 
long read alignments. 


Input data
----------

 * Ciliate MAC genome assembly, Fasta format.
 * PacBio HiFi or CCS read library mapping onto that assembly, sorted/indexed
   BAM format; mapper should report valid CIGAR string and NM tag.


Installation
------------

Dependencies are specified as a Conda environment YAML file `env.yaml`. Create a
Conda environment with the specified dependencies, then install bleties locally
with `pip`:

```bash
conda env create -f env.yaml -n bleties
conda activate bleties
cd /path/to/bleties # path to this folder
pip install -e .
```


Usage
-----

For help, use the `-h` or `--help` option, with or without the subworkflow 
names:

```
bleties --help
bleties milraa --help
bleties miser --help
bleties milret --help
bleties milcor --help
bleties miltel --help
bleties insert --help
```


In addition there are two scripts for plotting the output from the MILRAA and
MILCOR modules. See the `--help` messages for usage instructions:

```
milcor_plot.py --help
milraa_plot.py --help
```


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
 * [MILTEL](miltel.md) -- Identify potential chromosome breakage sites with
   telomere addition
 * Insert -- Utility for inserting/deleting IES sequences from a reference
   genome to generate MAC+IES or MAC-IES versions.

