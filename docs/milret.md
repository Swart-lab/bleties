MILRET - Method of IES Long-read Retention
==========================================

Inputs
------

 * Mapping of PacBio HiFi/CCS reads to reference genome
 * Reference assembly
 * Feature table (GFF3) of IES junctions, either from MILRAA output or third
   party tool

Reimplementation of MIRET in the ParTIES pipeline. MIRET takes the following
inputs: two reference genomes (somatic and germline), separate mapping files of
the same reads to the somatic and germline genomes, and coordinates of known
IES junctions. Reads that map with match to the somatic version are counted as
IES-, whereas reads that map with match to the germline version are counted as
IES+.

With long reads, we assume that IESs are spanned completely by most reads, vs
short reads, where reads are unlikely to completely span an IES insert. So we
do not count soft/hard clips on the ends of reads, only matches and inserts. We
also assume that the read mapper will handle mapping of reads containing inserts
properly. Therefore, we only map the reads to the somatic genome. Reads that do
contain IES+ forms will be reported as mappings with insert operations ("I" in
the CIGAR string). We then simply compare reads mapping with match to the
somatic genome at the IES junction (counted as IES-) to reads mapping with an
insert at the exact location of the IES junction (counted as IES+).

Differences of MILRET to MIRET:
 * MILRET uses only the somatic (IES-) as mapping reference, and assumes that 
   the mapper will correctly handle IES+ reads by mapping them with an insert 
   relative to the reference (I operation in the CIGAR string)
 * MILRET also counts and reports reads with additional deletions relative to
   the mapping reference. These may represent alternative excisions, 
   misassembly, or misalignments



