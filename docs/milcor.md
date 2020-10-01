MILCOR - Method of IES Long-read CORrelation
============================================

Inputs
------

 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed)
 * Feature table (GFF3) of IES junctions, either from MILRAA output or third
   party tool

With long reads (>1 kbp) it is possible to count IES retention at the level of
individual reads. In PacBio or Nanopore sequencing, libraries are prepared
without amplification (e.g. by PCR), so the reads represent original molecules,
without the possibility of PCR chimerism. We could therefore potentially
classify reads into MIC-origin or MAC-origin, in the case of vegetative cells,
or examine the dynamics of IES excision in developing MACs.

MILCOR reports a per-read IES retention score that complements the per-IES
retention score reported by MILRET. This is not possible with short read
sequencing where reads typically do not span an entire IES. In the calculation
of the per-IES retention score, reads that do not span at least one defined IES
junction site are not counted.
