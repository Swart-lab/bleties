#!/usr/bin/env python3

import pysam
import subprocess
import logging
from collections import defaultdict

from bleties.SharedFunctions import getCigarOpQuerySeqs, nested_dict_to_list_fixed_depth
from bleties.SharedFunctions import Gff


logger = logging.getLogger("Miltel")


def softclipped_seqs_from_bam(alnfile):
    """Get softclipped query sequence segments from BAM mappings

    Parameters
    ----------
    alnfile : pysam.AlignmentFile

    Returns
    -------
    list
        list of dicts, each dict representing a softclipped sequence segment,
        with the keys seq, qname, qlen, qstart, qend, rname, rstart, rend.
        Coordinates are zero-based end-exclusive (python convention)
    """
    out = []
    counter = 0
    for rec in alnfile.fetch():
        if (not rec.is_supplementary) and (not rec.is_unmapped) and (not rec.is_secondary):
            counter += 1
            if counter % 1000 == 0:
                logger.debug(f"Processed {str(counter)} entries")
            rname = rec.reference_name
            qname = rec.query_name
            qlen = rec.query_length # includes soft-clips but not hard-clips
            opseqs = getCigarOpQuerySeqs(
                        rec.query_sequence, # includes soft-clips but not hard-clips
                        rec.cigartuples,
                        rec.reference_start,
                        target_op="S")
            if opseqs and len(opseqs) > 0:
                for (seq, qstart, qend, rstart, rend) in opseqs:
                    out.append({
                        "seq": seq, 
                        "qname": qname,
                        "qlen": qlen,
                        "qstart": qstart, # coordinates 0-based [)
                        "qend": qend,
                        "rname": rname,
                        "rstart": rstart,
                        "rend": rend})
    return(out)


def find_telomeres(seq, telomere="ACACCCTA", minlength=24):
    """Find telomere sequences with NCRF in a list of sequences

    Assumes that NCRF is in the PATH.

    Parameters
    ----------
    seq : str
        Sequence to be scanned
    telomere : str
        Sequence of the telomere repeat. Default is ACACCCTA, corresponding to
        Blepharisma and Stentor telomere.
    minlength : int
        Minimum length of consecutive telomeres to detect (bp)
    """
    ncrf_out = subprocess.run(
            ["NCRF", f"telomere:{telomere}", f"--minlength={str(minlength)}"],
            capture_output=True,
            input=str.encode(seq))
    return(ncrf_out)


def parse_NCRF_output(ncrf):
    """Parse NCRF output and report alignment coordinates and orientation

    Assumes NCRF output without --positionalevents option

    Parameters
    ----------
    ncrf : byte
        Raw NCRF output

    Returns
    -------
    list
        List of dicts with the following fields: seqname, alnstart, alnend,
        orientation, score
    """
    if ncrf.stderr.decode().startswith("(0 alignments"):
        # No alignments in this segment
        return(None)
    else:
        lines = ncrf.stdout.decode().split("\n")
        # drop comment lines and split by whitespace
        lines = [line.split() for line in lines if not line.startswith("#")]
        lines = [line for line in lines if len(line) > 0]
        if len(lines) % 2 != 0:
            logger.error("Number of NCRF output lines is not a multiple of two")
        else:
            out = []
            for i in range(int(len(lines)/2)):
                seqname = lines[i*2][0]
                alncoords = lines[i*2][3] # e.g. 1234-2345
                orientation = lines[i*2 + 1][0][-1] # last character is + or -
                score = int(lines[i*2 + 1][2].split("=")[1]) # e.g. score=1234
                alnstart, alnend = [int(i) for i in alncoords.split("-")]
                out.append({"seqname": seqname,
                    "alnstart" : alnstart,
                    "alnend" : alnend,
                    "orientation" : orientation,
                    "score" : score})
            return(out)


def rekey_softclip_recs_by_ref(seqlist, telomere="ACACCCTA", min_telomere_length=24):
    """Rekey softclip_seqs_from_bam output to a dict organized by reference
    
    Search for telomeres in softclip seqs and parse NCRF output

    Parameters
    ----------
    seqlist : list
        list of dicts, output from softclip_seqs_from_bam()

    Returns
    -------
    dict
        dict keyed by rname -> rstart -> clip_orientation
        for soft clips, assume that rstart == rend
    """
    out = defaultdict( # rname
            lambda: defaultdict( # rstart
                lambda: defaultdict( # clip_orientation
                    list)))
    counter = 0
    for rec in seqlist:
        counter += 1
        if counter % 1000 == 0:
            logger.debug(f"Processed {str(counter)} entries")
        if rec['rstart'] == rec['rend']:
            orientation = None
            if rec['qstart'] == 0:
                orientation = "left"
            if rec['qend'] == rec['qlen']:
                orientation = "right"
            if orientation:
                # Find telomere sequence
                ncrf = find_telomeres(rec['seq'], telomere, min_telomere_length)
                ncrf_parse = parse_NCRF_output(ncrf)
                if ncrf_parse:
                    rec['ncrf_parse'] = ncrf_parse
                out[rec['rname']][rec['rstart']][orientation].append(rec)
            else:
                logger.warn(f"No clipping orientation found for softclip {str(rec)}")
        else:
            logger.warn(f"Softclip with nonzero extent on reference, {str(rec)}")

    return(out)


def call_features_from_rekeyed_dict(seqdict):
    """Call putative chromosome breakage sites junctions from rekeyed dict

    Parameters
    ----------
    seqdict : dict
        Output of rekey_softclip_recs_by_ref()

    Returns
    -------
    Gff
        Gff object reporting coordinates of chromosome breakage sites
    """
    out = Gff()
    reclist = nested_dict_to_list_fixed_depth(seqdict, 3)
    for rname, rstart, orientation, recs in reclist:
        out_aln_orientations = []
        out_aln_gaps = []
        for rec in recs:
            if "ncrf_parse" in rec:
                if orientation == "left":
                    # clip on left of read, so check distance from end of 
                    # alignment to the end of the softclipped segment
                    gap = len(rec['seq']) - rec['ncrf_parse'][0]['alnend'] # TODO case where there is more than one aligned segmnet
                elif orientation == "right":
                    # clip on right of read, so check distance from beginning
                    # of alignment to the start of softclipped segment
                    gap = rec['ncrf_parse'][0]['alnstart']
                out_aln_orientations.append(rec['ncrf_parse'][0]['orientation'])
                out_aln_gaps.append(gap)
        if len(out_aln_gaps) > 0:
            gffid = f"CBS_{rname}_{str(rstart+1)}_{orientation}"
            attrs=[f"ID={gffid}",
                   f"orientation={orientation}",
                   "telomere_senses=" + " ".join([str(i) for i in out_aln_orientations]),
                   "gaps_to_telomere=" + " ".join([str(i) for i in out_aln_gaps])]
            out.addEntry( # self, linearr, gffid
                [rname, "MILTEL", "chromeosome_breakage_site",
                 rstart+1, rstart+1, # Convert to 1-based coords for GFF
                 len(out_aln_gaps),
                 '.', '.',
                 ";".join(attrs)],
                 gffid)
    return(out)
