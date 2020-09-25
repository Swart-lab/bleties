To do list for BleTIES
======================

 - [ ] MILCOR - visualize correlations between IESs relative to distance
 - [ ] MILRAA - find suitable clustering and fuzzy matching parameters for
     PacBio subreads/CLR reads.
 - [ ] Implement Assembly module 
    - Will have to chop up the reads because PacBio reads are long, also want to
       include some of the context. 
    - Have to somehow deal with the potential mapping ambiguities at boundaries
       of inserts
 - [ ] Implement Insert/Remove module, to remove putative IESs from assembly, or
     add them back in
 - [ ] Implement internal telomere finder as separate module in BleTIES
