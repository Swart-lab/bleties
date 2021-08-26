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

If the genome has existing annotations in a GFF3 feature table, this can be
supplied to the `--featuregff` option. The coordinates of the annotations will
be updated after the insertion of IESs, and other annotation fields remain
unchanged. If a feature is interrupted by one or more IESs, the feature will be
split into two or more features. If the feature contains an ID attribute and
the option `--addsuffix` is chosen, then the ID of the original feature is
suffixed with `.seg_0`, `.seg_1` et seq. to distinguish the different segments
split from the original entry. By default the ID field is unchanged. The
updating of an existing feature table is only available in insert mode, because
the deletion or truncation of an existing feature annotation will require
manual intervention or curation.

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
