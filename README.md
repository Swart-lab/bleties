Blepharisma Toolbox for Interspersed DNA Elimination Studies (BleTIES)
======================================================================

![BLETIES logo](./docs/bleties_logo.png)

This is a reimplementation of [ParTIES](https://github.com/oarnaiz/ParTIES) for 
long read alignments. 

The required inputs are a ciliate MAC genome assembly, and a PacBio HiFi (high-
fidelity CCS reads) read library mapped onto that assembly. The mapper should
report a valid CIGAR string and NM tag (for number of mismatches) per aligned
read.


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

Scripts for plotting and visualizing data are in the `scripts/` subfolder.

Refer to the [full documentation](./docs/index.md) for more information.


Citations
---------

Please cite dependencies if you use them:

 * `pysam` - A Heger, K Jacobs, et al. [https://github.com/pysam-developers/pysam ](https://github.com/pysam-developers/pysam)
 * `htslib`, `samtools` - H Li, et al. 2009. ["The Sequence Alignment/Map format and SAMtools"](https://doi.org/10.1093/bioinformatics/btp352) _Bioinformatics_ 25 (16) : 2078-2079.
 * `biopython` - PJA Cock, et al. 2009. ["Biopython: freely available Python tools for computational molecular biology and bioinformatics"](https://doi.org/10.1093/bioinformatics/btp163) _Bioinformatics_ 25 (11) : 1422-1423.
 * `muscle` - RC Edgar, 2004. ["MUSCLE: multiple sequence alignment with high accuracy and high throughput"](https://doi.org/10.1093/nar/gkh340) _Nucleic Acids Research_ 32 (5) : 1792-1797.
 * `ncrf` - RS Harris, M Cechova, KD Makova, 2019. ["Noise-cancelling repeat finder: uncovering tandem repeats in error-prone long-read sequencing data"](https://doi.org/10.1093/bioinformatics/btz484) _Bioinformatics_ 35 (22) : 4809-4811.
 * `spoa` - R Vaser, I Sovic, N Nagarajan, Mile Sikic, 2017. ["Fast and accurate de novo genome assembly from long uncorrected reads"](https://doi.org/10.1101/gr.214270.116) _Genome Research_ 27 : 737-746.
