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
cd /path/to/bleties # path to this repo
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

BleTIES is research software. Please [cite us](CITATION.md) if you use the
software in a publication.
