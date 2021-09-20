Citing BleTIES and dependencies
===============================

BleTIES is research software. If you use it in a publication, please cite our
paper:

> Brandon K. B. Seah, Estienne C. Swart. (2021) BleTIES: Annotation of natural
> genome editing in ciliates using long read sequencing. *Bioinformatics*
> btab613; doi: https://doi.org/10.1093/bioinformatics/btab613

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
@article{10.1093/bioinformatics/btab613,
    author = {Seah, Brandon K B and Swart, Estienne C},
    title = "{BleTIES: Annotation of natural genome editing in ciliates using long read sequencing}",
    journal = {Bioinformatics},
    year = {2021},
    month = {09},
        abstract = "{Ciliates are single-celled eukaryotes that eliminate specific, interspersed DNA sequences (internally eliminated sequences, IESs) from their genomes during development. These are challenging to annotate and assemble because IES-containing sequences are typically much less abundant in the cell than those without, and IES sequences themselves often contain repetitive and low-complexity sequences. Long read sequencing technologies from Pacific Biosciences and Oxford Nanopore have the potential to reconstruct longer IESs than has been possible with short reads, but require a different assembly strategy. Here we present BleTIES, a software toolkit for detecting, assembling, and analyzing IESs using mapped long reads.BleTIES is implemented in Python 3. Source code is available at https://github.com/Swart-lab/bleties (MIT license), and also distributed via Bioconda.Benchmarking of BleTIES with published sequence data.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btab613},
    url = {https://doi.org/10.1093/bioinformatics/btab613},
    note = {btab613},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btab613/40316268/btab613.pdf},
    }
```
