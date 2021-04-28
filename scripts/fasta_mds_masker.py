#!/usr/bin/env python

#Example execution (first mask IESs, then mask MDSs with "N"):
#./gff_fasta_ies_masker.py mac_plus_ies.gff mac_plus_ies.fa > mac_plus_ies.ies_masked.fa
#./fasta_mds_masker.py mac_plus_ies.ies_masked.fa > mac_plus_ies.mds_masked.fa

from Bio import SeqIO, Seq
from sys import argv
import re

pat = re.compile("[ACGTN]")

for rec in SeqIO.parse(open(argv[1]), "fasta"):
  seq = pat.sub("N", str(rec.seq))
  print(f">{rec.name}\n{seq}")
    

# vim:sts=2:ts=2:sw=2
