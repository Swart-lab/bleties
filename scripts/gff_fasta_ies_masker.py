#!/usr/bin/env python

#Example execution:
#./gff_fasta_ies_masker.py mac_plus_ies.gff mac_plus_ies.fa > mac_plus_ies.ies_masked.fa

from Bio import SeqIO, Seq
from sys import argv

mask_d = {}
for line in open(argv[1]):
  if line[0] != '#':
    atoms = line.split("\t")
    mask_d.setdefault(atoms[0], []).append((int(atoms[3]) -1, int(atoms[4])))

for rec in SeqIO.parse(open(argv[2]), "fasta"):
  seq = str(rec.seq)
  try:
    ies_coords = mask_d[rec.name]
    for coord_pair in ies_coords:
      seq = seq[:coord_pair[0]] + seq[coord_pair[0]:coord_pair[1]].lower() + seq[coord_pair[1]:]
  
    print(f">{rec.name}\n{seq}")
    
  except KeyError:
    pass

# vim:sts=2:ts=2:sw=2
