#!/usr/bin/env python3

import pysam
import subprocess
import logging

from bleties.SharedFunctions import getCigarOpQuerySeqs


logger = logging.getLogger("Miltel")


def softclipped_seqs_from_bam(alnfile):
    """Get softclipped query sequence segments from BAM mappings

    Parameters
    ----------
    alnfile : pysam.AlignmentFile
    """
    out = []
    for rec in alnfile.fetch():
        if (not rec.is_supplementary) and (not rec.is_unmapped) and (not rec.is_secondary):
            rname = rec.reference_name
            qname = rec.query_name
            opseqs = getCigarOpQuerySeqs(
                        rec.query_sequence,
                        rec.cigartuples,
                        rec.reference_start,
                        target_op="S")
            if opseqs and len(opseqs) > 0:
                for (seq, qstart, qend, rstart, rend) in opseqs:
                    out.append((seq, 
                        qname, qstart, qend,
                        rname, rstart, rend))
    return(out)

def find_telomeres_in_softclipped_seqs(seqlist, telomere="ACACCCTA", minlength=24):
    """Find telomere sequences with NCRF in a list of sequences

    Assumes that NCRF is in the PATH.

    Parameters
    ----------
    seqlist : list
        Output from softclipped_seqs_from_bam(), a list of tuples, where the
        first element in each tuple is a str representing softclipped seq.
    telomere : str
        Sequence of the telomere repeat. Default is ACACCCTA, corresponding to
        Blepharisma and Stentor telomere.
    minlength : int
        Minimum length of consecutive telomeres to detect (bp)
    """
    out = []
    for (seq, qname, qstart, qend, rname, rstart, rend) in seqlist:
        ncrf_out = subprocess.run(
                ["NCRF", f"telomere:{telomere}", f"--minlength={str(minlength)}"],
                capture_output=True,
                input=str.encode(seq))
        out.append(ncrf_out)
    return(out)
