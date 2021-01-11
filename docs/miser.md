MISER - Method of IES Spurious or Erroneous Reporting
=====================================================

This is an experimental module, may fail in unexpected ways.

Inputs
------

 * Mapping of PacBio HiFi/CCS reads to reference genome (preferably BAM, sorted
   and indexed)
 * Reference assembly
 * Feature table of putative IES coordinates, produced by MILRAA

MISER takes an existing set of IES predictions, produced by MILRAA, and screens
it for potential mispredictions caused by paralogy, misassembly, or erroneous
mappings. 

For each putative IES (insertion or deletion), the set of reads mapping to that 
site is found, and split into two subsets: those containing the indel and those
without. For each subset, the mean percent mismatch of alignments vs. the 
reference is taken. 

 * If either subset has high (>5%) mismatch rate, "high error" is reported.
 * If the subset with indel has a significantly higher mismatch than the subset
   without, possible paralog is reported.
 * If the subset without indel has a significantly higher mismatch, possible 
   misassembly is reported.
 * Otherwise the putative IES is "ok".
