Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)
======================================================================

![BLETIES logo](./bleties_logo.png)

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for
long read alignments. 


Input data
----------

 * Ciliate MAC genome assembly, Fasta format.
 * PacBio read library mapping onto that assembly, sorted/indexed BAM format;
   mapper should report valid CIGAR string and NM tag.


Installation
------------

Dependencies are specified as a Conda environment YAML file `env.yaml`. Create a
Conda environment with the specified dependencies, then install bleties locally
with `pip`:

```bash
cd /path/to/bleties # path to this folder
conda env create -f env.yaml -n bleties
conda activate bleties
pip install .
```

Run tests after installation:

```bash
python -m unittest -v bleties.TestModule
```


Usage
-----

To list input arguments and options and their respective usage, use the `-h` or
`--help` option, with or without the subworkflow names:

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
   mapping of PacBio subreads or CCS reads to reference genome assembly.
 * [MISER](miser.md) -- Screen for potentially erroneous IES calls, to curate
   the list of putative IES junctions (experimental!).
 * [MILRET](milret.md) -- Use curated IES junctions, or IES junction
   coordinates from third-party tools, to calculate IES retention scores from
   mapping of CCS reads to a reference MAC assembly.
 * [MILCOR](milcor.md) -- Calculate per-read IES retention scores and bin reads
   as MIC/MAC in origin, from mapping of CCS reads to MAC assembly.
 * [MILTEL](miltel.md) -- Identify potential chromosome breakage sites with
   telomere addition, from mapping of CCS reads to MAC assembly.
 * [Insert](insert.md) -- Utility for inserting IES sequences from a MAC
   reference genome to generate MAC+IES hybrid sequence. Also can perform the
   reverse operation.


Citations
---------

Please cite dependencies if you use them:

 * `pysam` - A Heger, K Jacobs, et al. [https://github.com/pysam-developers/pysam ](https://github.com/pysam-developers/pysam)
 * `htslib`, `samtools` - H Li, et al. 2009. ["The Sequence Alignment/Map format and SAMtools"](https://doi.org/10.1093/bioinformatics/btp352) _Bioinformatics_ 25 (16) : 2078-2079.
 * `biopython` - PJA Cock, et al. 2009. ["Biopython: freely available Python tools for computational molecular biology and bioinformatics"](https://doi.org/10.1093/bioinformatics/btp163) _Bioinformatics_ 25 (11) : 1422-1423.
 * `muscle` - RC Edgar, 2004. ["MUSCLE: multiple sequence alignment with high accuracy and high throughput"](https://doi.org/10.1093/nar/gkh340) _Nucleic Acids Research_ 32 (5) : 1792-1797.
 * `ncrf` - RS Harris, M Cechova, KD Makova, 2019. ["Noise-cancelling repeat finder: uncovering tandem repeats in error-prone long-read sequencing data"](https://doi.org/10.1093/bioinformatics/btz484) _Bioinformatics_ 35 (22) : 4809-4811.
 * `spoa` - R Vaser, I Sovic, N Nagarajan, Mile Sikic, 2017. ["Fast and accurate de novo genome assembly from long uncorrected reads"](https://doi.org/10.1101/gr.214270.116) _Genome Research_ 27 : 737-746.
