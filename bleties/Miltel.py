#!/usr/bin/env python3

import pysam
import subprocess
import logging
import json
from collections import defaultdict

from bleties.SharedFunctions import getCigarOpQuerySeqs, nested_dict_to_list_fixed_depth
from bleties.SharedFunctions import report_summary_string, report_list_modes
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


class Miltel(object):
    """Object to store softclip accounting records"""

    def __init__(self, alnfile, refgenome):
        """Initialize Miltel object with alignment file and reference genome

        Parameters
        ----------
        alnfile : pysam.AlignmentFile
            Input BAM file opened in pysam
        refgenome : Bio.SeqRecords
            Reference genome file
        """
        self._alnfile = alnfile
        self._refgenome = refgenome
        self._clippedseqs = None # Dict to store clipped seq parses
        self._clippeddict = None # Dict to store clipped seq parses


    def dump(self):
        """Data dump of Miltel object in JSON format"""
        outstr = json.dumps({"clippedseqs": self._clippedseqs, 
                             "_clippeddict": self._clippeddict}, 
                            indent = 2)
        return(outstr)


    def get_softclips(self):
        """Parse BAM file for softclipped reads"""
        self._clippedseqs = softclipped_seqs_from_bam(self._alnfile)


    def find_telomeres(self, telomere="ACACCCTA", min_telomere_length=24):
        """Rekey softclip_seqs_from_bam output to a dict organized by reference
        
        Search for telomeres in softclip seqs and parse NCRF output

        This function should be called after get_softclips()

        Parameters
        ----------
        telomere : str
            Telomere sequence (5') to search
        min_telomere_length : int
            Minimum telomere alignment length to count

        Returns
        -------
        dict
            dict keyed by rname -> rstart -> clip_orientation
            for soft clips, assume that rstart == rend
        """
        self._clippeddict = defaultdict( # rname
                                lambda: defaultdict( # rstart
                                    lambda: defaultdict( # clip_orientation
                                        list)))
        counter = 0
        for rec in self._clippedseqs:
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
                    self._clippeddict[rec['rname']][rec['rstart']][orientation].append(rec)
                else:
                    logger.warn(f"No clipping orientation found for softclip {str(rec)}")
            else:
                logger.warn(f"Softclip with nonzero extent on reference, {str(rec)}")


    def report_CBS_GFF(self):
        """Call putative chromosome breakage sites junctions from rekeyed dict

        This function must be called after find_telomeres

        Returns
        -------
        Gff
            Gff object reporting coordinates of chromosome breakage sites
        """
        out = Gff()
        reclist = nested_dict_to_list_fixed_depth(self._clippeddict, 3)
        for rname, rstart, orientation, recs in reclist:
            out_aln_orientations = []
            out_aln_gaps = []
            for rec in recs:
                if "ncrf_parse" in rec:
                    gap = None
                    currgap = None
                    alnorientation = None
                    # If there is more than one alignment, find the one
                    # closest to the boundary
                    for alnrec in rec['ncrf_parse']:
                        if orientation == "left":
                            # clip on left of read, so check distance from end of 
                            # alignment to the end of the softclipped segment
                            currgap = len(rec['seq']) - alnrec['alnend']
                        elif orientation == "right":
                            # clip on right of read, so check distance from beginning
                            # of alignment to the start of softclipped segment
                            currgap = alnrec['alnstart']
                        if gap:
                            if currgap and currgap < gap:
                                gap = currgap
                                alnorientation = alnrec['orientation']
                        else:
                            gap = currgap
                            alnorientation = alnrec['orientation']
                    out_aln_orientations.append(alnorientation)
                    out_aln_gaps.append(gap)
            if len(out_aln_gaps) > 0:
                gffid = f"CBS_{rname}_{str(rstart+1)}_{orientation}"
                telomere_sense = "_".join(report_list_modes(out_aln_orientations))
                telomere_gap_average = round(sum(out_aln_gaps)/len(out_aln_gaps), 3)
                telomere_senses = report_summary_string(out_aln_orientations)
                telomere_gaps = report_summary_string(out_aln_gaps)
                attrs=[f"ID={gffid}",
                       f"orientation={orientation}",
                       f"telomere_sense={telomere_sense}",
                       f"telomere_gap_average={telomere_gap_average}",
                       f"telomere_senses={telomere_senses}",
                       f"telomere_gaps={telomere_gaps}"]
                out.addEntry( # self, linearr, gffid
                    [rname, "MILTEL", "chromosome_breakage_site",
                     rstart+1, rstart+1, # Convert to 1-based coords for GFF
                     len(out_aln_gaps),
                     '.', '.',
                     ";".join(attrs)],
                     gffid)
        return(out)
