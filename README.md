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

Refer to the [full documentation](./docs/index.md) for more information.
