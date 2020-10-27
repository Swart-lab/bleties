Insert
======

For most ciliate species, MAC genome data and reference assemblies are more
readily available than MIC genome data. For certain analyses, it is useful to
have a reference sequence containing IESs (IES+, or MAC+IES), if a MIC genome
is not available.

To facilitate this, the Insert utility module allows users to artificially
insert IES sequences into a MAC reference genome from a GFF3 file and a Fasta
of the corresponding IES sequences, to obtain the MAC+IES reference and a GFF3
file of the updated IES coordinates (`--mode insert`). The utility also permits
the reverse operation (`--mode delete`).


Input data
----------

The inputs required are a Fasta file containing the reference genome (supplied
to `--fasta` option), and GFF3 file of IES junction/region coordinates
(`--ies`).

In insert mode (`--mode insert`), an additional Fasta file of the IES sequences
is required, supplied to the `--iesfasta` option. The headers in this Fasta
file should correspond to the `ID` value in the attributes field of the IES GFF
file.


Output
------

If run in insert mode (`--mode insert`), IES sequences are inserted into the
MAC genome to produce a hybrid MAC+IES reference, with the following files:

 * `{OUT}.iesplus.fasta` - Fasta file containing MAC+IES reference
 * `{OUT}.iesplus.gff` - GFF file containing updated coordinates for the
   inserted IES sequences, derived from the input GFF file of IES junctions.

If run in delete mode (`--mode delete`), IES sequences are removed from the
MAC+IES (or MIC) genome to produce a MAC-like reference, with the following
files:

 * `{OUT}.iesminus.fasta` - Fasta file containing IES- reference
 * `{OUT}.iesminus.gff` - GFF file containing updated coordinates for the
   removed IES sequences, derived from the input GFF file of IES regions.
