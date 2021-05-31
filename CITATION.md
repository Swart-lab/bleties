Citing BleTIES and dependencies
===============================

BleTIES is research software. If you use it in a publication, please cite our
preprint:

> Brandon K. B. Seah, Estienne C. Swart. (2021, preprint) BleTIES: Annotation
> of natural genome editing in ciliates using long read sequencing. *bioRxiv*
> 2021.05.18.444610; doi: https://doi.org/10.1101/2021.05.18.444610 

The software repository itself can be cited with either the GitHub URL or the
Zenodo archive DOI:

[![DOI](https://zenodo.org/badge/294123134.svg)](https://zenodo.org/badge/latestdoi/294123134)

In addition, you should also cite the following dependencies that BleTIES is
built upon: 

 * `pysam` - A Heger, K Jacobs, et al. [https://github.com/pysam-developers/pysam ](https://github.com/pysam-developers/pysam)
 * `htslib`, `samtools` - H Li, et al. 2009. ["The Sequence Alignment/Map format and SAMtools"](https://doi.org/10.1093/bioinformatics/btp352) _Bioinformatics_ 25 (16) : 2078-2079.
 * `biopython` - PJA Cock, et al. 2009. ["Biopython: freely available Python tools for computational molecular biology and bioinformatics"](https://doi.org/10.1093/bioinformatics/btp163) _Bioinformatics_ 25 (11) : 1422-1423.
 * `muscle` - RC Edgar, 2004. ["MUSCLE: multiple sequence alignment with high accuracy and high throughput"](https://doi.org/10.1093/nar/gkh340) _Nucleic Acids Research_ 32 (5) : 1792-1797.
 * `ncrf` - RS Harris, M Cechova, KD Makova, 2019. ["Noise-cancelling repeat finder: uncovering tandem repeats in error-prone long-read sequencing data"](https://doi.org/10.1093/bioinformatics/btz484) _Bioinformatics_ 35 (22) : 4809-4811.
 * `spoa` - R Vaser, I Sovic, N Nagarajan, Mile Sikic, 2017. ["Fast and accurate de novo genome assembly from long uncorrected reads"](https://doi.org/10.1101/gr.214270.116) _Genome Research_ 27 : 737-746.


BibTeX
------

```
@article {Seah2021.05.18.444610,
	author = {Seah, Brandon K. B. and Swart, Estienne C.},
	title = {BleTIES: Annotation of natural genome editing in ciliates using long read sequencing},
	elocation-id = {2021.05.18.444610},
	year = {2021},
	doi = {10.1101/2021.05.18.444610},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Summary Ciliates are single-celled eukaryotes that eliminate specific, interspersed DNA sequences (internally eliminated sequences, IESs) from their genomes during development. These are challenging to annotate and assemble because IES-containing sequences are much less abundant in the cell than those without, and IES sequences themselves often contain repetitive and low-complexity sequences. Long read sequencing technologies from Pacific Biosciences and Oxford Nanopore have the potential to reconstruct longer IESs than has been possible with short reads, and also the ability to detect correlations of neighboring element elimination. Here we present BleTIES, a software toolkit for detecting, assembling, and analyzing IESs using mapped long reads.Availability and implementation BleTIES is implemented in Python 3. Source code is available at https://github.com/Swart-lab/bleties (MIT license), and also distributed via Bioconda.Contact Contact: kb.seah{at}tuebingen.mpg.deSupplementary information Benchmarking of BleTIES with published sequence data.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2021/05/18/2021.05.18.444610},
	eprint = {https://www.biorxiv.org/content/early/2021/05/18/2021.05.18.444610.full.pdf},
	journal = {bioRxiv}
	}
```
