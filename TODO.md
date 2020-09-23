To do list for MILRAA
=====================

 - [ ] Look for Ds that directly abut I inserts, as proxy for alternative
    excision boundaries?
 - [ ] Question: How many putative IES junctions are spanned by each read?
    - How "chimeric" is a given read? I.e. if a read contains non-excised IESs,
       does it have more non-excised IESs on average?
    - Per-read IES retention score - a measure of "linkage" between IES
       retention events
 - [x] Group putative IESs that are overlapping into clusters, to identify
     potential alternative excisions
 - [x] Cluster putative IESs so that they do not have to exactly overlap / have
     exactly the same length
   - [x] Cluster IES inserts by sequence distance
 - [x] Adjust putative IES junctions ...
   - [x] Adjust coordinates so that the junctions are TA junctions, if possible
   - [x] Adjust coordinates to maximize length of pointer, if possible
 - [ ] Implement Assembly module 
    - Will have to chop up the reads because PacBio reads are long, also want to
       include some of the context. 
    - Have to somehow deal with the potential mapping ambiguities at boundaries
       of inserts
 - [ ] Implement MIRET (MILRET) module to report IES retention score
    - [x] For each junction, report M operations spanning that junction 
        -- Representing IES- form
    - [x] Also report D operations spanning that junction. If there are many D
        operations spanning that junction, something odd is going on!
    - [x] Report I operations at that junction, of all insert lengths
        -- Representing IES+ form
    - [ ] ? Report soft/hard clips that touch that junction (for completeness'
        sake)
    - [ ] Look more carefully into how MIRET implements the IES+ and IES-
        counting, and why they find it necessary to use mapping to both germline
        and somatic genomes
    - Should provide a separate shortlist of IES breakpoints, which can be from
        MILRAA, but should not directly get the result from MILRAA. This is
        because user may already have a curated set of IES junctions and is just
        seeking a retention score calculation, vs. de novo prediction of IESs.
    - Should also report how many potential IESs are not in the provided
        shortlist. (Actually this is better done by just running MILRAA separately
        and comparing the results.)
    - ParTIES MIRET uses mapping to both MAC and MIC genomes to compare counts
        of reads spanning IES- version vs. reads spanning IES+ version. With
        PacBio reads, and known IES junctions, why don't we just rely on the
        Insert operation of the mapper? 
 - [ ] Implement Insert/Remove module, to remove putative IESs from assembly, or
     add them back in
 - [ ] MILRAA prediction of inserts: how to deal with fuzzy boundaries and
     variable insert lengths? Will be more of a problem with subreads vs. CCS
     reads. Maybe give user option to specify "error profile" of input reads
 - [x] Unit tests
 - [ ] Implement internal telomere finder as separate module in BleTIES
 - [x] Sanity check - check that contigs in reference match headers in BAM
