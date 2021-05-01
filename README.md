Basic Long-read Enabled Toolkit for Interspersed DNA Elimination Studies (BleTIES)
==================================================================================

![BLETIES logo](./docs/bleties_logo.png)

[![DOI](https://zenodo.org/badge/294123134.svg)](https://zenodo.org/badge/latestdoi/294123134)
![Bioconda](https://img.shields.io/conda/vn/bioconda/bleties.svg)
![License](https://img.shields.io/github/license/Swart-lab/bleties.svg)

BleTIES is a tool for prediction and targeted assembly of internally eliminated
sequences (IESs) in ciliate genomes, using single-molecule long read
sequencing. The design and name of the software was inspired by
[ParTIES](https://github.com/oarnaiz/ParTIES).

The required inputs are a ciliate MAC genome assembly, and a long read
sequencing library (PacBio subreads or error-corrected CCS reads, or Nanopore
reads) mapped onto that assembly. The mapper should report a valid CIGAR string
and NM tag (for number of mismatches) per aligned read.


Installation
------------


### Install released version with Conda

The released versions are distributed via Bioconda, and can be installed with Conda:

```bash
# Create new environment called "bleties"
conda create -c conda-forge -c bioconda -n bleties bleties
# Activate environment
conda activate bleties
# Check version and view help message
bleties --version
bleties --help
# Run tests
python -m unittest -v bleties.TestModule
```


### Install development version

If you want to test the latest development version, clone this Git repository,
then install with pip.

Dependencies are specified as a Conda environment YAML file `env.yaml`. Create a
Conda environment with the specified dependencies, then install with `pip`:

```bash
git clone git@github.com:Swart-lab/bleties.git
cd bleties
conda env create -f env.yaml -n bleties_dev
conda activate bleties_dev
pip install .
```

Run tests after installation:

```bash
python -m unittest -v bleties.TestModule
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

Two scripts are available for plotting visualizations of MILRAA and MILCOR
results: `milcor_plot.py` and `milraa_plot.py`.

Additional useful scripts are in the folder `scripts/`

Refer to the [full documentation](./docs/index.md) for more information.


Citations
---------

BleTIES is research software. Please [cite us](CITATION.md) if you use the
software in a publication.
