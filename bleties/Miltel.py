#!/usr/bin/env python3

import pysam
import subprocess
import logging
import json
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from bleties.SharedFunctions import getCigarOpQuerySeqs, nested_dict_to_list_fixed_depth
from bleties.SharedFunctions import report_summary_string, report_list_modes
from bleties.SharedFunctions import Gff
from bleties.Milraa import alnFromSeqs


logger = logging.getLogger("Miltel")


def softclipped_seqs_from_bam(alnfile, ctg=None, start=None, stop=None):
    """Get softclipped query sequence segments from BAM mappings

    Parameters
    ----------
    alnfile : pysam.AlignmentFile
        Alignment to be processed
    ctg : str
        Name of contig to be processed. If None, then process entire alignment.
    start : int
    stop : int
        Coordinates in the contig to be processed. If None, then process entire
        contig. 0-based python coordinates.

    Returns
    -------
    list
        list of dicts, each dict representing a softclipped sequence segment,
        with the keys seq, qname, qlen, qstart, qend, rname, rstart, rend.
        Coordinates are zero-based end-exclusive (python convention)
    """
    out = []
    counter = 0
    for rec in alnfile.fetch(ctg, start, stop):
        if (not rec.is_supplementary) and (not rec.is_unmapped) and (not rec.is_secondary):
            counter += 1
            if counter % 1000 == 0:
                logger.debug(f"Processed {str(counter)} entries")
            rname = rec.reference_name
            qname = rec.query_name
            qlen = rec.query_length  # includes soft-clips but not hard-clips
            opseqs = getCigarOpQuerySeqs(
                rec.query_sequence,  # includes soft-clips but not hard-clips
                rec.cigartuples,
                rec.reference_start,
                target_op="S")
            if opseqs and len(opseqs) > 0:
                for (seq, qstart, qend, rstart, rend) in opseqs:
                    out.append({
                        "seq": seq,
                        "qname": qname,
                        "qlen": qlen,
                        "qstart": qstart,  # coordinates 0-based [)
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
            logger.error(
                "Number of NCRF output lines is not a multiple of two")
        else:
            out = []
            for i in range(int(len(lines)/2)):
                seqname = lines[i*2][0]
                alncoords = lines[i*2][3]  # e.g. 1234-2345
                orientation = lines[i*2 + 1][0][-1]  # last character is + or -
                score = int(lines[i*2 + 1][2].split("=")[1])  # e.g. score=1234
                alnstart, alnend = [int(i) for i in alncoords.split("-")]
                out.append({"seqname": seqname,
                            "alnstart": alnstart,
                            "alnend": alnend,
                            "orientation": orientation,
                            "score": score})
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
        self._clippedseqs = None  # Dict to store clipped seq parses
        self._clippeddict = None  # Dict to store clipped seq parses

    def dump(self):
        """Data dump of Miltel object in JSON format"""
        outstr = json.dumps({"clippedseqs": self._clippedseqs,
                             "_clippeddict": self._clippeddict},
                            indent=2)
        return(outstr)

    def get_softclips(self, ctg=None, start=None, stop=None):
        """Parse BAM file for softclipped reads"""
        # TODO: do this for only a defined region of the alignment
        self._clippedseqs = softclipped_seqs_from_bam(
            self._alnfile, ctg, start, stop)

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
        self._clippeddict = defaultdict(  # rname
            lambda: defaultdict(  # rstart
                lambda: defaultdict(  # clip_orientation
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
                    ncrf = find_telomeres(
                        rec['seq'], telomere, min_telomere_length)
                    ncrf_parse = parse_NCRF_output(ncrf)
                    if ncrf_parse:
                        rec['ncrf_parse'] = ncrf_parse
                    self._clippeddict[rec['rname']
                                      ][rec['rstart']][orientation].append(rec)
                else:
                    logger.warn(
                        f"No clipping orientation found for softclip {str(rec)}")
            else:
                logger.warn(
                    f"Softclip with nonzero extent on reference, {str(rec)}")

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
                telomere_sense = "_".join(
                    report_list_modes(out_aln_orientations))
                telomere_gap_average = round(
                    sum(out_aln_gaps)/len(out_aln_gaps), 3)
                telomere_senses = report_summary_string(out_aln_orientations)
                telomere_gaps = report_summary_string(out_aln_gaps)
                readcov = self._alnfile.count(
                    str(rname), start=int(rstart), stop=int(rstart)+1,
                    read_callback=lambda r: (not r.is_supplementary) and (not r.is_secondary))
                breakscore = None
                if readcov > 0:
                    # Round to 3 sig figs
                    breakscore = float('%.3g' % float(len(out_aln_gaps)/readcov))
                # attr.append(f"average_coverage={str(readcov)}")
                attrs = [f"ID={gffid}",
                         f"orientation={orientation}",
                         f"telomere_sense={telomere_sense}",
                         f"telomere_gap_average={telomere_gap_average}",
                         f"telomere_senses={telomere_senses}",
                         f"telomere_gaps={telomere_gaps}",
                         f"average_coverage={str(readcov)}"]
                out.addEntry(  # self, linearr, gffid
                    [rname, "MILTEL", "chromosome_breakage_site",
                     rstart+1, rstart+1,  # Convert to 1-based coords for GFF # TODO verify that +1 should not apply here because softclipping is not ref consuming
                     str(breakscore),
                     '.', '.',
                     ";".join(attrs)],
                    gffid)
        return(out)

    def report_other_clips_GFF_fasta(self, min_clip_length=50):
        """Call other clipped sequences and call consensus of their sequences

        This function must be called after find_telomeres

        Parameters
        ----------
        min_clip_length : int
            Minimum length of clipped sequence to consider

        Returns
        -------
        Gff
            Gff object reporting coordinates of clip junctions
        list
            List of SeqRecord objects reporting clipped sequences (consensus
            sequences if more than one sequence at a site). SeqRecord IDs
            correspond to ID fields in Gff output
        """
        out_gff = Gff()
        out_seqs = []
        reclist = nested_dict_to_list_fixed_depth(self._clippeddict, 3)
        for rname, rstart, orientation, recs in reclist:
            gffid = f"CLIPJUNCTION_{rname}_{str(rstart+1)}_{orientation}"
            seqs = []
            counter = 0
            for rec in recs:
                # Skip telomere sequences, and those that are too short
                if (not "ncrf_parse" in rec) and (len(rec['seq']) > min_clip_length):
                    seqs.append(SeqRecord(Seq(rec['seq']),
                                          id=f"{gffid}_{str(counter)}"))
                    counter += 1
            if len(seqs) == 1:
                out_seqs.extend(seqs)
            elif len(seqs) > 1:
                # Align with Muscle and get consensus sequence
                # TODO decide if better to have dumb or gap consensus for our purposes
                logger.debug(
                    f"Consensus from {str(len(seqs))} seqs for {gffid}")
                seqs_cons = alnFromSeqs(seqs)
                seqs_cons.id = f"{gffid}_cons"
                out_seqs.append(seqs_cons)

            if len(seqs) > 0:
                attrs = [f"ID={gffid}",
                         f"orientation={orientation}"]
                out_gff.addEntry(  # self, linearr, gffid
                    [rname, "MILTEL", "clip_junction",
                     rstart+1, rstart+1,  # Convert to 1-based coords for GFF # TODO verify that +1 should not apply here because softclipping is not ref consuming
                     len(seqs),  # Number of clipped segments over threshold
                     '.', '.',
                     ";".join(attrs)],
                    gffid)

        return(out_gff, out_seqs)
