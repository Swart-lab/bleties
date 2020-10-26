#!/usr/bin/env python3

import re
import json
from collections import defaultdict
from statistics import median
import pysam
import logging
import subprocess
from io import StringIO
from tempfile import NamedTemporaryFile
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio import AlignIO

from bleties.SharedValues import SharedValues
from bleties.SharedFunctions import *


# Define logger
logger = logging.getLogger("Milraa")


# TODO: refactor with SharedFunctions.getCigarOpQuerySeqs
def getIndels(cigar, pos, minlength, qseq):
    """Parse cigar string and alignment position and report insertions or
    deletions.

    Return list of tuples (start pos, end pos, insert length,
    insertion or deletion, insert sequence). If it is a deletion then no insert
    length or insert sequence are reported because that can be parsed from the
    reference.

    Parameters
    ----------
    cigar : str
        CIGAR string of alignment to process
    pos : int
        Alignment start position on reference, 1-based numbering
    minlength : int
        Minimum length of indel to report
    qseq : str
        Query sequence

    Returns
    -------
    list
        list of tuples (int, int, int, str, str) representing (start pos, end
        pos, insert length, insertion/deletion, insert sequence)
    """
    outarr = []  # Array to store tuples of results
    # Initialize counters for how much query or reference seq is consumed
    ref_consumed = 0
    que_consumed = 0
    # Split cigar string into individual operations
    cigs = re.findall(r"\d+[\w\=]", cigar)
    for cig in cigs:
        cigmatch = re.match(r"(\d+)([\w\=])", cig)  # Get number and operation
        # We want to look for insert operations above a min length,
        # map their positions on the reference, and count how many reads
        # support a given putative insert
        # This requires that we count the operations that consume reference
        # and add it to the POS field
        if cigmatch.group(2) == "I":  # If insert operation,
            ins_len = int(cigmatch.group(1))  # Length of the current insert
            if ins_len >= minlength:  # Check that insert is above min length
                # Get the start and end positions of the insert
                # 1-based, insert is to the right of position
                ins_pos_start = pos + ref_consumed - 1
                # Get the sequence of the insert
                # 0-based, following pysam convention
                ins_seq = qseq[que_consumed:que_consumed + ins_len]
                outarr.append(
                    (ins_pos_start, ins_pos_start, ins_len, "I", ins_seq))
        # We also look for delete operations above a min length
        # These are already present in the reference, so the "insert length" is 0
        if cigmatch.group(2) == "D":  # If delete operation,
            del_len = int(cigmatch.group(1))  # Length of current deletion
            if del_len >= minlength:
                # Get start and end pos of deletion
                del_pos_start = pos + ref_consumed  # 1-based, deletion starts ON this position
                del_pos_end = pos + ref_consumed + del_len - \
                    1  # 1-based, deletion ends ON this position
                # If deletion, no insert sequence reported
                outarr.append((del_pos_start, del_pos_end, 0, "D", ""))
        # Count ref and query consumed _after_ the insert has been accounted for
        if cigmatch.group(2) in SharedValues.REFCONSUMING:
            ref_consumed += int(cigmatch.group(1))
        if cigmatch.group(2) in SharedValues.QUERYCONSUMING:
            que_consumed += int(cigmatch.group(1))
    return(outarr)


def getIndelJunctionSeqs(iesgff, iesconsseq, ref, flanklen):
    """Get sequence at indel junctions.

    Parameters
    ----------
    iesgff : Gff
        Gff object listing insertions and deletions (output from
        reportPutativeIes). Coordinates in Gff are 1-based inclusive.
    iesconsseq : dict
        dict of SeqRecords for indels (output from reportPutativeIes)
    ref : dict
        Reference genome (dict of SeqRecord objects)
    flanklen : int
        Length of flanking sequence at junctions to report

    Returns
    -------
    list
        list of list (breakpoint, left junction seq, right junction seq)
    """
    outseqs = []
    for breakpointid in iesgff:
        ctg = iesgff.getValue(breakpointid, 'seqid')
        start = int(iesgff.getValue(breakpointid, 'start'))
        end = int(iesgff.getValue(breakpointid, 'end'))
        flankleftseq = ""
        flankrightseq = ""
        indel = ""
        refunedited = ""
        # Check if this is insertion junction or deletion region
        # GFF convention is to record zero-length features with start=end, and
        # junction site is to right of coordinate.
        if start == end:  # insertion junction
            indel = "I"
        elif start < end:
            indel = "D"
        else:
            raise Exception(
                f"Start cannot be less than end, feature {breakpointid}")
        # Get the left flanking junction on reference
        # Check whether sequence with flanking will run off the end
        if start >= flanklen:  # TODO Check off-by-one
            if indel == "I":  # For insert, feature starts on RIGHT of coordinate
                flankleftseq = ref[ctg][start - flanklen: start].seq.lower()
            elif indel == "D":  # For deletion, feature starts ON the coordinate
                flankleftseq = ref[ctg][start-1-flanklen:start-1].seq.lower()
        else:  # Pad the left side
            flankleftseq = (flanklen - start) * "-" + \
                ref[ctg][0:start].seq.lower()
        # Get the right flanking junction on reference
        # Check whether sequence with flanking will run off the end
        if (end + flanklen) <= len(ref[ctg]):
            flankrightseq = ref[ctg][end: end + flanklen].seq.lower()
        else:
            flankrightseq = ref[ctg][end:].seq.lower()
        # Add the IES sequence fragment
        if indel == "I":
            flankleftseq = flankleftseq + \
                iesconsseq[breakpointid][0:flanklen].seq.upper()
            flankrightseq = iesconsseq[breakpointid][-flanklen:].seq.upper() + \
                flankrightseq
        elif indel == "D":
            # In GFF, start and end are 1-based and BOTH inclusive
            flankleftseq = flankleftseq + \
                ref[ctg][start-1: start-1+flanklen].seq.upper()
            flankrightseq = ref[ctg][end -
                                     flanklen:end].seq.upper() + flankrightseq
        # Get the reference sequence
        if indel == "I":
            refunedited = ref[ctg][start-flanklen:end+flanklen].seq.lower()
        elif indel == "D":
            # Start minus one because for deletion, start is ON the deleted region
            refunedited = ref[ctg][start-flanklen-1:end+flanklen].seq.lower()
        # Find pointers if present
        (pointer, pointerstart, pointerend) = getPointers(
            ref[ctg], start, end, iesconsseq[breakpointid], breakpointid)
        if pointerstart != start:
            logger.info(
                f"Position of pointer {breakpointid} has been adjusted")
            start = pointerstart
            end = pointerend
        # Convert pointers to TA junctions if possible
        (tastart, taend, tapointer) = adjustPointerTA(
            pointerstart, pointerend, pointer)
        # Maximize pointer lengths
        (ppstart, ppend, pppointer) = adjustPointerMaxlength(
            ref[ctg], start, end, pointer, iesconsseq[breakpointid])
        # Put everything together and return
        outseqs.append([breakpointid,
                        str(ctg),
                        str(start),
                        str(end),
                        str(flankleftseq),
                        str(flankrightseq),
                        str(pointer),
                        str(tastart),
                        str(taend),
                        str(tapointer),
                        str(ppstart),
                        str(ppend),
                        str(pppointer),
                        str(iesconsseq[breakpointid].seq),
                        str(refunedited)
                        ])
    return(outseqs)


def getPointers(seq, start, end, iesseq, name):
    """Find potential pointer sequences at putative IES junctions

    Report the pointer sequences as-is from the mapper/GFF file

    Parameters
    ----------
    seq : str
        Sequence of contig with putative IES junction/region or SeqRecord
    start : int
        Start position of putative IES junction/region, 1-based (from GFF)
    end : int
        End position of putative IES junction/region, 1-based (from GFF). If
        start == end , the IES is not present on the reference sequence and the
        junction is to the RIGHT of position, per GFF convention.
        If start < end, the IES is retained in the reference sequence.
    iesseq : str or SeqRecord
        In case where start == end, IES sequence must be supplied separately
    name : str
        ID or name of breakpoint, for error reporting.

    Returns
    -------
    str
        Putative pointer sequence
    int, int
        Adjusted IES start and IES end positions, such that the pointer sequence
        in the MDS is to the _right_ of the insert.
    """
    # Give provisional name to this indel if none is supplied
    if not name:
        name = f"{seq} {int(start)} {int(end)}"

    if start == end:
        indel = "I"
    elif start < end:
        indel = "D"
        iesseq = seq[start-1: end]
    else:
        raise Exception(f"Start cannot be less than end, feature {name}")

    # Remove gap characters from ies sequence
    if isinstance(iesseq, str):
        iesseq = iesseq.replace("-", "")
    elif isinstance(iesseq, SeqRecord):
        iesseq = iesseq.seq.ungap("-")
    else:
        raise Exception(
            f"iesseq must be of type str or SeqRecord but is {type(iesseq)}")

    pointer = ""
    pointerstart, pointerend = start, end
    leftcheck = ""
    rightcheck = ""
    # check left of IES
    i = 0
    while iesseq[i] == seq[end+i]:
        leftcheck += iesseq[i]
        i += 1
        # break out of loop if the next position runs off edge
        if i >= len(iesseq) or end+i >= len(seq):
            break
    # check right of IES
    i = -1
    if indel == "I":  # beacuse GFF puts zero-length features to right of coordinate
        refstart = start
    elif indel == "D":
        refstart = start - 1
    while iesseq[i] == seq[refstart+i]:
        rightcheck = iesseq[i] + rightcheck
        i -= 1
        # break out of loop if next position runs off edge
        if -i > len(iesseq) or refstart+i < 0:
            break
    # take the longer putative pointer
    if len(leftcheck) < 2 and len(rightcheck) < 2:
        pointer = None
    elif len(leftcheck) > len(rightcheck) and len(leftcheck) >= 2:
        pointer = leftcheck
    elif len(rightcheck) > len(leftcheck) and len(rightcheck) >= 2:
        pointer = rightcheck
        # If the insert position is to the right of the putative pointer seq,
        # report adjusted pointer position such that the insert position always
        # lies to the left of the putative pointer seq on the MDS.
        pointerstart = pointerstart - len(pointer)
        pointerend = pointerend - len(pointer)
    elif len(rightcheck) == len(leftcheck):
        logger.info(f"Breakpoint {name} has potential pointers on both sides")
        pointer = "tie"  # if both sides could potentially have a pointer
    else:
        logger.warn(
            f"Unexpected result in pointer search for breakpoint {name}")
    return(pointer, pointerstart, pointerend)


def adjustPointerTA(start, end, pointer):
    """Adjust putative pointer to TA junction

    If a putative pointer sequence contains TA, it could be a TA junction.
    Adjust junction coordinates to make it a TA junction if possbile.
    If there is more than one TA in the pointer, take the first occurrence.

    Parameters
    ----------
    start : int
    end : int
        Same as getPointers()
    pointer : seq
        Pointer sequence, reported from getPointers()

    Returns
    -------
    int, int
        Adjusted start and end coordinates, such that TAs lie to the right of
        the junction coordinates. If no TA is present, None is returned
    str
        Adjusted pointer sequence beginning with TA. If no TA is present, None
        is returned.
    """
    tastart, taend = None, None
    pointernew = None
    if pointer:
        mm = re.search(r"TA", pointer)
        if mm:
            displace = mm.start()
            pointernew = pointer[displace:]
            tastart = start + displace
            taend = end + displace
    return(tastart, taend, pointernew)


def adjustPointerMaxlength(seq, start, end, pointer, iesseq):
    """Maximize pointer length by checking upstream to try and extend the
    repeat sequence

    The read mapper may already be reporting indel junction coordinates that
    have been optimized in this way, but we do this check because different
    mappers may behave differently. The pointer sequence should lie to the
    right of the junction coordinate on the MDS.

    Parameters
    ----------
    seq : str
        Contig reference sequence
    start : int
    end : int
    pointer : str
        Same as for adjustPointerTA()
    iesseq : str
        IES sequence, necesary if start == end

    Returns
    -------
    int, int
        Adjust start and end coordinates of junction
    str
        Adjusted pointer sequence
    """
    # Get IES sequence for a deletion feature
    if start < end:
        indel = "D"
        iesseq = seq[start - 1: end]
    elif start == end:
        indel = "I"
        if not iesseq:
            raise Exception("IES sequence is missing")
    # Look upstream of the pointer to try and extend match
    i = -1
    if indel == "I":  # beacuse GFF puts zero-length features to right of coordinate
        refstart = start
    elif indel == "D":
        refstart = start - 1
    if not pointer:
        pointer = ""
    while iesseq[i] == seq[refstart+i]:
        pointer = iesseq[i] + pointer
        i -= 1
        # break out of loop if next position runs off edge
        if -i > len(iesseq) or refstart+i < 0:
            break
    # Return adjusted coordinates
    return(start + 1 + i, end + 1 + i, pointer)


def alignSeqsMuscle(seqlist, muscle_path="muscle"):  # TODO: Use SPOA instead?
    """Align list of sequences with Muscle and return alignment

    Parameters
    ----------
    seqlist : list
        List of SeqRecord objects
    muscle_path : str
        Path to Muscle binary

    Returns
    -------
    MultipleSeqAlignment
        Alignment of the input sequences
    """
    from Bio.Align.Applications import MuscleCommandline
    from io import StringIO
    hh = StringIO()
    SeqIO.write(seqlist, hh, "fasta")
    data = hh.getvalue()
    muscle_cl = MuscleCommandline(muscle_path, clwstrict=True)
    aln_stdout, aln_stderr = muscle_cl(stdin=data)
    aln = AlignIO.read(StringIO(aln_stdout), "clustal")
    return(aln)


def alnFromSeqs(seqlist, threshold=0.7):
    """Align list of sequences with Muscle and report consensus

    Parameters
    ----------
    seqlist : list
        list of str representing sequences, may be different lengths
    threshold : float
        Threshold for consensus, passed to gap_consensus()

    Returns
    -------
    SeqRecord
        Consensus sequence, including gaps
    """
    if isinstance(seqlist[0], str):
        seqrecs = [SeqRecord(Seq(i, generic_dna)) for i in seqlist]
    elif isinstance(seqlist[0], Seq):
        seqrecs = [SeqRecord(i) for i in seqlist]
    elif isinstance(seqlist[0], SeqRecord):
        seqrecs = seqlist
    else:
        raise Exception(
            "sequence list must comprise str, Seq, or SeqRecord objects")
    aln = alignSeqsMuscle(seqrecs)
    alninf = AlignInfo.SummaryInfo(aln)
    alncons = alninf.gap_consensus(threshold=threshold)
    alnconsrec = SeqRecord(alncons)
    return(alnconsrec)


def alnDumbFromSeqs(seqlist, threshold=0.7):
    """Report consensus of a list of sequences that are all the same length

    Parameters
    ----------
    seqlist : list
        list of str, representing sequences all of the same length
    threshold : float
        Threshold for consensus, passed to dumb_consensus()

    Returns
    -------
    SeqRecord
        "Dumb" consensus sequence
    """
    if isinstance(seqlist[0], str):
        # If not a SecRecord object, transform to SeqRecords
        seqrecs = [SeqRecord(Seq(i, generic_dna)) for i in seqlist]
    elif isinstance(seqlist[0], Bio.Seq.Seq):
        seqrecs = [SeqRecord(i) for i in seqlist]
    elif isinstance(seqlist[0], Bio.SeqRecord.SeqRecord):
        seqrecs = seqlist
    else:
        raise Exception(
            "sequence list must comprise str, Seq, or SeqRecord objects")
    aln = MultipleSeqAlignment(seqrecs)
    alninf = AlignInfo.SummaryInfo(aln)
    alncons = alninf.dumb_consensus(threshold=threshold)
    alnconsrec = SeqRecord(alncons)
    return(alnconsrec)


def spoaConsensus(seqs, mode=1):
    """Get consensus sequence from a set of PacBio subreads (or subread
    fragments) with SPOA

    Parameters
    ----------
    seqs : list
        list of SeqRecord objects, representing PacBio subreads or subread
        fragments
    mode : int
        Algorithm mode for SPOA, recommended 1 (global). Passed to `-l`
        argument.

    Returns
    -------
    SeqRecord
        Consensus sequence
    """
    fh = NamedTemporaryFile(suffix=".fasta", mode="w", delete=False)
    SeqIO.write(seqs, fh.name, "fasta")
    fh.close()
    cons = subprocess.run(
        ["spoa", f"-l {mode}", "-r 0", fh.name], capture_output=True).stdout.decode()
    cons = SeqIO.read(StringIO(cons), "fasta")
    os.remove(fh.name)
    return(cons)


def findLongestInsert(seqs, rname, rstart=0):
    """Find insert corresponding to longest gap in the reference sequence

    Parameters
    ----------
    seqs : dict
        Dict of two aligned SeqRecord objects, the insert-containing called 
        "Consensus" and the reference called rname
    rstart : int
        Original start coordinate of the reference sequence

    Returns
    -------
    str
        Insert sequence
    int
        Coordinate of the inserted read, with respect to the original rstart,
        0-based, the insert lies to the RIGHT of the coordinate
    """
    gaphits = [i for i in re.finditer(r"\-+", str(seqs[rname].seq))]
    longestgap_coords = sorted(gaphits, reverse=True, key=lambda x: x.span()[
                               1] - x.span()[0])[0].span()
    # Find gaps in ref sequence before the longest insert
    offset = sum([i.end() - i.start()
                  for i in gaphits if i.end() < longestgap_coords[0]])
    # subtract 1 because insert should be to the right of the coordinate
    insert_start = rstart + longestgap_coords[0] - offset - 1
    return(str(seqs['Consensus'].seq[longestgap_coords[0]:longestgap_coords[1]]),
           insert_start)

# TODO: complain if the 2nd largest gap is more than 50% the length of the
# longest gap


def subreadNamesToZmwCoverage(qnames):
    """From list of PacBio subread names, report number of ZMWs represented

    QNAME of a PacBio subread has the following convention:
    {movieName}/{holeNumber}/{qStart}_{qEnd}
    We want to count the number of holes (ZMWs), because a given hole may result
    in multiple subreads.

    Parameters
    ----------
    qnames : list
        read names of PacBio subreads

    Returns
    -------
    int
        Number of ZMWs represented by the above subreads
    """
    zmwnames = ["/".join(n.split("/")[0:2]) for n in qnames]
    zmwnames = set(zmwnames)
    return(len(zmwnames))


def alnRemoveGapOnlyCols(aln):
    """Remove gap-only columns from a Bio.Align object

    Parameters
    ----------
    aln : Bio.Align
        Alignment to check

    Returns
    -------
    Bio.Align
        Alignment without the gap only columns
    """
    alninfo = AlignInfo.SummaryInfo(aln)
    alncons = str(alninfo.gap_consensus(threshold=1))
    gapcols = [hit.span() for hit in re.finditer(r"\-+", alncons)]
    notgaps = get_not_gaps(0, len(alncons), gapcols)
    alnnogaps = aln[:,notgaps[0][0]:notgaps[0][1]]
    if len(notgaps) > 1:
        for i in range(1,len(notgaps)):
            alnnogaps = alnnogaps + aln[:, notgaps[i][0]:notgaps[i][1]]
    return(alnnogaps)


class IesRecords(object):
    """Records of putative IESs from mappings"""

    def __init__(self, alnfile, alnformat, refgenome):
        """Constructor for IesRecords

        Internally represented by:
        _insSeqDict -- dict of sequences of detected inserts/deletions. Keys:
            contig (str) -> startpos (int) -> endpos (int) -> indel len (int) ->
            list of dicts containing data on the mapped read segments and their
            coordinates. The coordinates startpos and endpos in the keys are 1-
            based GFF convention. The coordinates in the dicts are 0-based
            python convention
        _alnfile -- as below
        _alnformat -- as below
        _refgenome -- as below

        Parameters
        ----------
        alnfile : pysam.AlignmentFile
            Alignment to parse
        alnformat : str
            Format of the alignment, either "bam" or "sam"
        refgenome : dict
            Reference genome sequences, dict of SeqRecord objects keyed by ID
        """
        # dict to store sequences of detected inserts/deletions
        self._insSeqDict = defaultdict(      # contig
            lambda: defaultdict(         # startpos
                lambda: defaultdict(     # endpos
                    lambda: defaultdict(  # indel length
                        list)            # list of dicts
                )
            )
        )
        # pysam.AlignmentFile object representing the BAM mapping
        self._alnfile = alnfile
        # Format of the alignment "bam" or "sam"
        self._alnformat = alnformat
        # Genome sequence used as reference for the mapping
        self._refgenome = refgenome

    def __str__(self):
        """Report summary stats of IesRecords object"""
        insseqdictlen = len(self._insSeqDict)
        nref = self._alnfile.nreferences
        alnformat = self._alnformat
        mapped = self._alnfile.mapped
        return("bleties.IesRecords object with "
               + "insSeqDict of length "
               + str(insseqdictlen)
               + " and alignment of format "
               + str(alnformat)
               + " with "
               + str(nref)
               + " references and "
               + str(mapped)
               + " mapped reads")

    def dump(self):
        """Data dump of IesRecords._insSeqDict in JSON format"""
        outstr = json.dumps(
            {"_insSeqDict": self._insSeqDict}, sort_keys=True, indent=2)
        return(outstr)

    def _addIndelsFromCigar(self, alignedsegment, minlength):
        """Check if alignment contains indels above minimum length, and record
        the corresponding breakpoints relative to the reference, and the insert
        length. 

        If the indel is an insert, insert length > 0. 
        If the indel is a deletion, insert length = 0. 

        Recorded in self._insSeqDict, keyed by contig -> start pos -> end pos
        -> insert length -> dict, where dict has fields
        'indelseq','qstart','qend','qname','rstart','rend','rname'.

        The coordinates in the dict keys are 1-based GFF type, whereas the
        coordinates in the dict are 0-based python type.

        Parameters
        ----------
        alignedsegment : pysam.AlignedSegment
            Record of aligned read from BAM file
        minlength : int
            Minimum length of indel for it to be recorded
        """
        # Look for inserts that are completely spanned by the read (i.e. I operations)
        rname = alignedsegment.reference_name
        qname = alignedsegment.query_name
        ins_tuples = getCigarOpQuerySeqs(alignedsegment.query_sequence,
                                         alignedsegment.cigartuples,
                                         alignedsegment.reference_start,
                                         "I")
        del_tuples = getCigarOpQuerySeqs(alignedsegment.query_sequence,
                                         alignedsegment.cigartuples,
                                         alignedsegment.reference_start,
                                         "D")
        if ins_tuples:
            for (indelseq, qstart, qend, rstart, rend) in ins_tuples:
                if qend - qstart >= minlength:
                    # Convert to 1-based end-inclusive for GFF
                    # indelstart=indelend for inserts
                    # it is correct to not +1 below
                    indelstart = rstart
                    indelend = rend
                    indellen = len(indelseq)
                    record = {'indelseq': indelseq,
                              'qstart': qstart, 'qend': qend, 'qname': qname,
                              'rstart': rstart, 'rend': rend, 'rname': rname}  # using 0-based coordinates
                    self._insSeqDict[rname][indelstart][indelend][indellen].append(
                        record)

        if del_tuples:
            for (indelseq, qstart, qend, rstart, rend) in del_tuples:
                if rend - rstart >= minlength:
                    # convert coords to 1-based end-inclusive numbering for GFF
                    # python uses 0-based end-exclusive, so no +1 for rend
                    indelstart = rstart + 1
                    indelend = rend
                    indellen = rend - rstart
                    record = {'indelseq': indelseq,  # this will be a blank sequence
                              'qstart': qstart, 'qend': qend, 'qname': qname,
                              'rstart': rstart, 'rend': rend, 'rname': rname}  # using 0-based coordinates
                    self._insSeqDict[rname][indelstart][indelend][indellen].append(
                        record)

    def findPutativeIes(self, minlength, ctg=None, start=None, stop=None):
        """Search alignment for clips and indels to identify putative IESs.
        Record them in the _insSeqDict

        Parameters
        ----------
        minlength : int
            Record only putative IESs of this length and above
        ctg : str
            Name of contig in reference to fetch alignments from
        start : int
        stop : int
            Coordinates in contig to fetch alignments from, 1-based, GFF-style
        """
        # convert coordinates to pysam 0-based
        if start:
            start = start - 1
        for line in self._alnfile.fetch(contig=ctg, start=start, stop=stop):
            if (not line.is_unmapped) and (not line.is_secondary) and (not line.is_supplementary):
                # total_mismatch = line.get_tag("NM") # Get number of mismatches # TODO record mismatches?
                # Find indels (putative IESs) over the minimum length and record them
                self._addIndelsFromCigar(line, minlength)

    def reportPutativeIesInsertFuzzy(self, mininsbreaks, mindelbreaks, dist_threshold=0.05):
        """Report putative IESs where insert lengths do not match exactly

        Parameters
        ----------
        mininsbreaks
        mindelbreaks
            Same as reportPutativeIes()
        dist_threshold : float
            Distance threshold for clustering, in range (0,1)

        Returns
        -------
        SharedFunctions.Gff
            Putative IES report parsable as GFF format
        dict
            dict of consensus sequences for putative IESs
        """
        # Initialize output
        gff = Gff()
        outseq = {}
        counter = 0

        # Cluster inserts (junctions)
        for rec in nested_dict_to_list_fixed_depth(self._insSeqDict, 3):
            [ctg, ins_start, ins_end, dd] = rec
            if ctg not in self._refgenome:
                logger.error(f"Sequence {ctg} not found in reference genome")

            breakpointid = None

            if ins_end == ins_start:
                ins_seqs = []
                for i in dd:
                    # Gather all sequences to be clustered
                    ins_seqs.extend([rec['indelseq'] for rec in dd[i]])
                clust_seqs, clust_ids = get_clusters_from_seqlist(
                    ins_seqs, dist_threshold)
                # For each cluster, report a putative IES
                for clust in clust_seqs:
                    counts = defaultdict(int)
                    for i in clust:
                        counts[len(i)] += 1
                    totalcount = sum(counts.values())
                    maxcounts = [l for l in counts if counts[l]
                                 == max(counts.values())]

                    # Report only if above minimum
                    if totalcount >= mininsbreaks:
                        # Put together ID for this breakpoint
                        if len(counts.keys()) == 1:  # cluster of one length
                            prefix = f"BREAK_POINTS"
                        elif len(counts.keys()) > 1:
                            prefix = f"BREAK_POINTS_FUZZY"
                        # Add median length to the bpid
                        breakpointfields = [
                            prefix, ctg, ins_start, ins_end, 
                            median([l for l in counts.keys()])]
                        breakpointid = "_".join([str(i)
                                                 for i in breakpointfields])
                        gfftype = "internal_eliminated_sequence_junction"

                        # Attributes list of key-value pairs
                        # report modal IES length
                        attr = [f"ID={breakpointid}",
                                "IES_length="+"_".join([str(l) for l in maxcounts])]
                        # Report number of counts per insert length
                        attr.append(
                            "cigar=" + " ".join([str(l) + "I*" + str(counts[l]) for l in counts]))

                        # Get indel consensus
                        if len(counts.keys()) == 1:  # cluster of one length
                            # take dumb consensus since all seqs are same length
                            consseq = alnDumbFromSeqs(clust)
                        else:
                            # Otherwise use Muscle (potentially slower) to align
                            consseq = alnFromSeqs(clust)

            elif ins_end > ins_start:
                del_len = ins_end - ins_start + 1  # GFF is end-inclusive
                totalcount = len(
                    self._insSeqDict[ctg][ins_start][ins_end][del_len])
                if totalcount >= mindelbreaks:
                    breakpointid = "_".join(["BREAK_POINTS", str(
                        ctg), str(ins_start), str(ins_end), str(del_len)])
                    gfftype = "internal_eliminated_sequence"
                    # get sequence from ref genome
                    indelseq = str(
                        self._refgenome[ctg].seq[ins_start - 1:ins_end])
                    consseq = SeqRecord(Seq(indelseq, generic_dna))
                    # Build attributes field
                    attr = [f"ID={breakpointid}",
                            f"IES_length={str(del_len)}"]
                    attr.append(f"cigar={str(del_len)}D*{str(totalcount)}")

            else:
                raise Exception("feature cannot start after it ends")

            # Record to Gff
            if breakpointid:
                # Get average coverage of region of interest
                if self._alnformat == "bam":
                    readcov = self._alnfile.count(
                        str(ctg), start=int(ins_start) - 1, stop=int(ins_end),
                        read_callback=lambda r: (not r.is_supplementary) and (not r.is_secondary))
                    attr.append(f"average_coverage={str(readcov)}")

                # Provisional approximate IES retention score
                # R = IES+ / (IES+ + IES-)
                # here the denominator is simply average coverage
                if readcov:
                    if gfftype == "internal_eliminated_sequence_junction":
                        provscore = round(totalcount/readcov, 4)
                    elif gfftype == "internal_eliminated_sequence":
                        # here deletions are counted
                        provscore = round((readcov-totalcount)/readcov, 4)

                # Find pointers if present
                ins_start, ins_end, pointerdict = self.reportAdjustPointers(
                    ctg, ins_start, ins_end, consseq, breakpointid)
                for i in pointerdict:
                    attr.append(f"{i}={str(pointerdict[i])}")

                consseq.id = breakpointid
                consseq.description = ";".join(attr)+";"
                outseq[consseq.id] = consseq

                outarr = [str(ctg),
                          "MILRAA",
                          gfftype,
                          str(ins_start),
                          str(ins_end),
                          str(provscore),
                          ".",
                          ".",
                          ";".join(attr)+";"]
                gff.addEntry(outarr, None)
                counter += 1
                if counter % 1000 == 0:
                    logger.info(f"Processed {counter} entries...")

        return(gff, outseq)
        # TODO Fuzzy cluster both the indel positions on ref and the ins lengths
        # Deletions: fuzzy cluster start and end positions separately.
        # This is somewhat trickier, because of the way the data are structured

    def reportAdjustPointers(self, ctg, ins_start, ins_end, consseq, breakpointid):
        """For a given indel sequence, find pointers and check for TA junctions

        Adjust to maximize pointer length if possible

        Parameters
        ----------
        ctg : str
            Name of sequence in self._refgenome
        ins_start : int
            Start position of indel, 1-based (GFF convention)
        ins_end : int
            End position of indel, 1-based inclusive (GFF convention). If
            equal to ins_start, then this is an insert junction.
        consseq : str or SeqRecord
            Consensus sequence of indel. Must be supplied if this is an insert.
        breakpointid : str
            Name for current indel

        Returns
        -------
        int, int
            ins_start and ins_end, which may be different from the input values
            if the pointer position has been adjusted
        dict
            Dict of pointer sequences and coordinates, also for putative TA
            junctions and maximized pointers, if present
        """
        out = {}
        # Find pointers if present
        (pointer, pointerstart, pointerend) = getPointers(
            self._refgenome[ctg], ins_start, ins_end, consseq, breakpointid)
        if pointerstart != ins_start:
            logger.info(
                f"Position of pointer {breakpointid} has been adjusted")
            ins_start = pointerstart
            ins_end = pointerend
        if pointer:  # Add pointer seq to attributes field if present
            out['pointer_seq'] = pointer
        # Convert pointers to TA junctions if possible
        (tastart, taend, tapointer) = adjustPointerTA(
            pointerstart, pointerend, pointer)
        # Add TA pointer sequence to attributes field if present
        if tapointer:
            out['ta_pointer_seq'] = tapointer
            out['ta_pointer_start'] = tastart
            out['ta_pointer_end'] = taend
        # Maximize pointer lengths, report if different from original
        (ppstart, ppend, pppointer) = adjustPointerMaxlength(
            self._refgenome[ctg], ins_start, ins_end, pointer, consseq)
        if pppointer and pppointer != pointer:
            out['pp_pointer_seq'] = pppointer
            out['pp_pointer_start'] = ppstart
            out['pp_pointer_end'] = ppend
        return(ins_start, ins_end, out)

    def reportPutativeIes(self, mininsbreaks, mindelbreaks):
        """After clips and indels have been recorded, report putative IESs above
        the minimum coverage, and if the input alignment is a BAM file, then
        also report average coverage in the breakpoint region. Output is written
        to open file handle. Min coverage for deletion breakpoints is expected
        to be lower because we are mapping to somatic genome in the typical use
        case, and reads with alternative excisions are thought to be rare.

        Parameters
        ----------
        mininsbreaks : int
            Minimum breakpoint coverage to report potential insertion
        mindelbreaks : int
            Minimum breakpoint coverage to report potential deletion

        Returns
        -------
        SharedFunctions.Gff
            Putative IES report parsable as GFF format
        dict
            dict of consensus sequences for putative IESs
        """
        # Create lists to hold SeqRecord objects and GFF output
        gff = Gff()
        outseq = {}

        # Parse the dict and report putative IESs above min coverage
        # We only check breakpoints which are completely spanned by a read ("I" or "D" operations)
        # TODO Reduce code duplication here
        for rec in nested_dict_to_list_fixed_depth(self._insSeqDict, 4):
            [ctg, ins_start, ins_end, ins_len, dd] = rec
            if ctg not in self._refgenome:
                logger.error(f"Sequence {ctg} not found in reference genome")
            evidencetype = None
            if ins_start == ins_end:
                evidencetype = "I"
            else:
                evidencetype = "D"
            countvalue = len(dd)
            breakpointid = None
            attr = []
            indel_len = 0

            # If the breakpoint is an insert type
            if evidencetype == "I" and countvalue >= mininsbreaks:
                breakpointid = "_".join(["BREAK_POINTS", str(
                    ctg), str(ins_start), str(ins_end), str(ins_len)])
                indel_len = ins_len
                gfftype = "internal_eliminated_sequence_junction"
                # Prepare attributes list of key-value pairs
                attr = [f"ID={breakpointid}",
                        f"IES_length={str(ins_len)}"]
                attr.append(f"cigar={str(ins_len)}I*{str(countvalue)}")

            # If the breakpoint is a deletion type
            elif evidencetype == "D" and countvalue >= mindelbreaks:
                # Add 1 because both start and end are inclusive
                del_len = int(ins_end) - int(ins_start) + 1
                indel_len = del_len
                breakpointid = "_".join(["BREAK_POINTS", str(ctg),
                    str(ins_start), str(ins_end), str(del_len)])
                gfftype = "internal_eliminated_sequence"
                # Build attributes field
                attr = [f"ID={breakpointid}",
                        f"IES_length={str(del_len)}"]
                attr.append(f"cigar={str(del_len)}D*{str(countvalue)}")

            # If the breakpoint has been defined, report it to the GFF file
            if breakpointid and attr:
                # Get read coverage from BAM file; SAM does not allow random access
                if self._alnformat == "bam":
                    readcov = self._alnfile.count(
                        str(ctg), start=int(ins_start) - 1, stop=int(ins_end),
                        read_callback=lambda r: (not r.is_supplementary) and (not r.is_secondary))
                    attr.append(f"average_coverage={str(readcov)}")

                # Provisional approximate IES retention score
                # R = IES+ / (IES+ + IES-)
                # here the denominator is simply average coverage
                if readcov:
                    if gfftype == "internal_eliminated_sequence_junction":
                        provscore = round(countvalue/readcov, 4)
                    elif gfftype == "internal_eliminated_sequence":
                        provscore = round((readcov-countvalue)/readcov, 4)

                # Get indel consensus
                consseq = self.reportIndelConsensusSeq(
                    ctg, ins_start, ins_end, indel_len)
                consseq.id = breakpointid
                consseq.description = ";".join(attr)+";"
                outseq[consseq.id] = consseq

                # Find pointers if present
                ins_start, ins_end, pointerdict = self.reportAdjustPointers(
                    ctg, ins_start, ins_end, consseq, breakpointid)
                for i in pointerdict:
                    attr.append(f"{i}={str(pointerdict[i])}")

                # Build GFF entry
                outarr = [str(ctg),            # 1 seqid
                          "MILRAA",            # 2 source
                          gfftype,             # 3 type
                          str(ins_start),      # 4 start
                          str(ins_end),        # 5 end
                          # 6 score - provisional IES retention score
                          str(provscore),
                          ".",                 # 7 strand
                          ".",                 # 8 phase
                          ";".join(attr)+";"   # 9 attributes
                          ]
                gff.addEntry(outarr, None)

        return(gff, outseq)

    def reportIndelConsensusSeq(self, ctg, indelstart, indelend, indellen):
        """Report consensus of indel sequence

        Parameters
        ----------
        ctg : str
            Name of contig
        indelstart : int
            Start position, 1-based
        indelend : int
            End position, 1-based
        indellen : int
            Length of indel

        Returns
        -------
        Bio.SeqRecord
            Consensus sequence of given indel
        """
        if indelstart == indelend:
            # From list of sequences as str, make a list of SeqRecord objects
            seqrecs = [SeqRecord(Seq(i['indelseq'], generic_dna))
                       for i in self._insSeqDict[ctg][indelstart][indelend][indellen]]
            # Make a pseudo-alignment from the list of SeqRecord objects
            aln = MultipleSeqAlignment(seqrecs)
            # Summarize alignment, and get dumb consensus sequence
            alninf = AlignInfo.SummaryInfo(aln)
            alncons = alninf.dumb_consensus()
            # alnconsrec = SeqRecord(alncons, id=consname)
            alnconsrec = SeqRecord(alncons)
            return(alnconsrec)
        elif int(indelend) > int(indelstart):
            # Get sequence from reference gneome
            indelseq = str(self._refgenome[ctg].seq[indelstart - 1:indelend])
            return(SeqRecord(Seq(indelseq, generic_dna)))

    def reportIndelConsensusSeqFuzzy(self, ctg, indelstart, indelend, indellens):
        """Report consensus alignment of insert sequence when length of insert
        is fuzzy matched

        Parameters
        ----------
        ctg : str
            Name of contig
        indelstart : int
            Start position, 1-based
        indelend : int
            End position, 1-based
        indellens : list
            list of ints corresponding to length of indel
        """

        # Gather all the sequences for all the lengths
        seqrecs = []
        for indellen in indellens:
            seqrecs.extend(
                [SeqRecord(Seq(i['indelseq'], generic_dna))
                 for i in self._insSeqDict[ctg][indelstart][indelend][indellen]])
        # Use Muscle to align these sequences
        aln = alignSeqsMuscle(seqrecs)
        # Summarize alignment, and get dumb consensus sequence
        alninf = AlignInfo.SummaryInfo(aln)
        # Take "dumb" consensus, but allowing gaps. Default is 70% majority consensus
        # We allow gaps because the dumb_consensus() function will take majority
        # base at a column as the consensus, even if most sequences at that
        # position have a gap
        alncons = alninf.gap_consensus()
        alnconsrec = SeqRecord(alncons)
        return(alnconsrec)

    def reportIndelReadMismatchPc(
            self, ctg, indelstart, indelend, indellen, minlength):
        """Report sequence mismatch % of query reads containing indel at a
        specific position vs. reads without indel.
        This is to flag indels that may originate from paralogs and hence are
        probably not true IESs.

        Parameters
        ----------
        ctg : str
            Name of contig
        indelstart : int
            Start position of indel, 1-based
        indelend : int
            End position, 1 based
        indellen : int
            Length of indel
        minlength : int
            Minimum length of putative IES insert to allow

        Returns
        -------
        list, list
            Two lists:
            ins_mm -- list of floats, mismatch % of query reads containing indel
                    at target position
            non_mm -- list of floats, mismatch % of query reads without indel at pos
        """
        # GFF allows start==end, but pysam does not recognise
        dummyend = indelend
        if indelend == indelstart:
            dummyend += 1
        # Get segments that overlap indel of interest
        # minus 1 for pysam uses 0-based coords
        itrr = self._alnfile.fetch(ctg, indelstart - 1, dummyend - 1)
        # Segments aligning to position of interest
        segs = [seg for seg in itrr]
        # initialize lists to report mismatch percentages
        non_mm = []
        ins_mm = []
        # for each segment, check if it contains indel at position of interest
        for seg in segs:
            indels = getIndels(seg.cigarstring, int(
                seg.reference_start), minlength, seg.query_sequence)
            indelcoords = set([(int(indel[0]), int(indel[1]))
                               for indel in indels])
            # number of mismatchs / query length * 100 pc
            mismatch_pc = 100 * float(seg.get_tag("NM")) / \
                float(seg.query_length)
            if (indelstart - 1, indelend - 1) in indelcoords:
                ins_mm.append(mismatch_pc)
            else:
                non_mm.append(mismatch_pc)
        return(ins_mm, non_mm)

    def getJuncClustersSeqs(self, clusters, rname, minseqs=10):
        """Get reference positions and alignment details for a set of
        coordinate clusters

        Parameters
        ----------
        clusters : list
            Output of SharedFunctions.get_clusters(), representing clusters of
            insert junction coordinates
        rname : str
            Name of reference contig
        minseqs : int
            Minimum number of subreads in a cluster to report. If this is too
            low, then the consensus insert sequence will not be accurate.

        Returns
        -------
        list
            List of dicts with keys 'positions' and 'seqs', corresponding
            values are lists containing insert coordinates (int) and dicts from
            the insSeqDict that contain info on the insert sequence.
        """
        jcs = []
        min_seqs = 10  # Minimum number of subreads in cluster to read
        for i in clusters:
            seqs = []
            for j in i:
                for jseqs in self._insSeqDict[rname][j][j].values():
                    seqs.extend(jseqs)
            if len(seqs) >= min_seqs:
                jcs.append({"positions": i, "seqs": seqs})
        return(jcs)

    def extractFlankingFromJcs(self, jcs, rname,
                               margin=100, lower_threshold=None, upper_threshold=None):
        """Extract insert of interest and flanking region from mapped reads

        Parameters
        ----------
        jcs : dict
            Element from the list output from getJuncClustersSeqs() above
        rname : str
            Name of contig of interest
        margin : int
            Sequence length (bp) to left and right of the insert to extract, to
            help anchor the consensus sequence against the reference
        lower_threshold : float
            Minimum insert length, expressed as fraction of the median insert
            length. Sequences below this are possibly fragments or sequencing
            or mapping errors and are ignored.
        upper_threshold : float
            Maximum insert length, expressed as multiple of the median insert
            length. Sequences above this may be sequencing or mapping errors
            and are ignored.

        Returns
        -------
        list
            list of SeqRecord of the extracted insert sequences and the flanking
            regions from mapped query sequences, Seq.id is the query QNAME
        tuple
            range of coordinates (min, max) on reference contig where the
            extracted insert sequences are positioned
        """
        coords = (min(jcs['positions']) - 1, max(jcs['positions']))
        details = jcs['seqs']
        alns = self._alnfile.fetch(rname, coords[0], coords[1])

        if lower_threshold and upper_threshold:
            # Keeping only sequnences within a defined range of median insert length
            lenmedian = median([len(i['indelseq']) for i in jcs['seqs']])
            lowthreshd = lower_threshold * lenmedian
            uppthreshd = upper_threshold * lenmedian
            details = [i for i in jcs['seqs'] if len(
                i['indelseq']) < uppthreshd and len(i['indelseq']) > lowthreshd]
            if len(details) == 0:
                logger.info(
                    f"Length thresholding discarded all inserts at {rname} {str(coords[0])}-{str(coords[1])}")
                logger.info("... insert lengths could be bimodally distributed")

        if len(details) > 0:
            # Key by query name
            subsetqnames = {i['qname']: {
                'qstart': i['qstart'], 'qend': i['qend']} for i in details}
            coords = (min([i['rstart'] for i in details]),
                      max([i['rend'] for i in details]))
            extractedseqs = []
            extractedids = []
            for i in alns:
                if (not i.is_supplementary) and (not i.is_secondary) and (not i.is_unmapped):
                    if i.query_name in subsetqnames:
                        start = max(
                            [subsetqnames[i.query_name]['qstart'] - margin, 0])
                        end = min([subsetqnames[i.query_name]
                                   ['qend'] + margin, i.query_length])
                        insertplusmargins = i.seq[start: end]
                        if len(insertplusmargins) == 0:
                            logger.debug(f"problem at {i.query_name}")
                            logger.debug(i.cigarstring)
                        extractedseqs.append(insertplusmargins)
                        extractedids.append(i.query_name)
            extractedrecs = [SeqRecord(
                Seq(extractedseqs[i]), id=extractedids[i]) for i in range(len(extractedseqs))]
            return(extractedrecs, coords)
        else:
            return(None, None)

    def spoaConsensusToFlanking(self, seqs, rname, rstart, rend, margin=100, mode=1):
        """Get consensus of subreads and align to flanking regions on reference
        genome

        Use SPOA to generate consensus sequence of subreads (or subread
        fragments) containing an insert relative to the reference, and align
        the consensus against flanking regions on reference genome contig, in
        order to identify the exact coordinate of the insert segment.

        If the consensus with insert is shorter than the original reference
        sequence without insert, then return None because it is likely to be an
        artefact.

        Parameters
        ----------
        seqs : list
            List of SeqRecord objects representing subread fragments that map to 
            region of interest and contain insert.
        rname : str
            Name of the contig of interest
        rstart : int
        rend : int
            Start and end coordinates on reference contig
        margin : int
            Length (bp) to left and right of rstart and rend to extract from
            reference sequence to align against the subreads consensus.

        Returns
        -------
        dict
            dict of SeqRecord objects keyed by SeqRecord.id field. IDs should
            be rname (given above) representing the reference sequence extract,
            and "Consensus" representing consensus of subread fragments
        """
        cons = spoaConsensus(seqs, mode=mode)
        # account for inserts close to ends of contig
        start = max([rstart - margin, 0])
        end = min([rend + margin, len(self._refgenome[rname].seq)])
        # Check that consensus of reads with insert is longer than original
        # without insert
        if len(cons) > (end - start):
            extract = self._refgenome[rname].seq[start:end]
            seqobj = SeqRecord(extract, id=rname)
            fh = NamedTemporaryFile(suffix=".fasta", mode="w", delete=False)
            SeqIO.write([seqobj, cons], fh.name, "fasta")
            fh.close()
            # key differences in settings: gap opening and extension penalties have
            # been changed to encourage long continuous gaps, which are what we
            # want to find when aligning IES+ consensus vs. IES- reference
            cons2 = subprocess.run(["spoa",
                                    "-e 0",  # gap extension penalty (default -6)
                                    "-g -16",  # gap opening penalty (default -8)
                                    "-r 1",  # alignment mode (global)
                                    fh.name],
                                   capture_output=True).stdout.decode()
            cons2 = SeqIO.parse(StringIO(cons2), "fasta")
            cons2 = MultipleSeqAlignment(cons2)
            # SPOA v4.0.3 leaves trailing gap-only cols in aln; remove these
            cons2 = alnRemoveGapOnlyCols(cons2)
            cons2 = {i.id: i for i in cons2}  # each i is a SeqRecord object
            os.remove(fh.name)
            return(cons2)
        else:
            logger.info(
                f"Inserts at {rname} {str(rstart)}-{str(rend)} anomalous; flanking + insert shorter than ref")
            return(None)

    def reportPutativeIesInsertSubreads(self, mininsbreaks, mindelbreaks,
                                        margin=100, max_cluster_dist=2,
                                        len_threshold=0.25):
        # TODO implement for deletions
        """Report putative IES inserts from subreads

        This function is called after self._insSeqDict has been produced by the
        function findPutativeIes

        Deletions are not yet implemented

        Parameters
        ----------
        mininsbreaks
        mindelbreaks
            Same as reportPutativeIes()
        margin : int
            Length of flanking sequence (bp) to extract around insert position
            from reads for assembling consensus sequence, to anchor the
            consensus afterwards by alignment to the reference.
        max_cluster_dist : int
            The first step is to look for inserts reported by the mapper that
            cluster together. This parameter sets the maximum distance (in bp)
            between insert coordinates to consider as part of the same cluster.
        len_threshold : float
            For the extracted sequences used to generate consensus, keep only
            those that are within a specific range of the median length. This
            is similar to the length filtering in PacBio CCS pipeline. The
            value should be a float between 0 and 1. For example if
            len_threshold is 0.25, then only sequences between 0.75 and 1.25
            times the median are used to generate consensus.

        Returns
        -------
        SharedFunctions.Gff
            Putative IES report parsable as GFF format
        dict
            dict of consensus sequences for putative IESs
        """
        # Initialize output
        gff = Gff()
        gff_problem = Gff()
        outseq = {}
        outseq_problem = {}
        counter = 0

        # For each contig, cluster putative IES inserts and get consensus seqs
        for rname in self._insSeqDict:
            # Find insert junctions (start == end) from insSeqDict
            # Find clusters that are no more than 2 bp apart
            juncs = [pos for pos in self._insSeqDict[rname]
                     if pos in self._insSeqDict[rname][pos]]
            juncs = [int(i) for i in juncs]
            # + 1 to max_cluster_dist because get_clusters is end-exclusive
            clusters = get_clusters(juncs, "bp", max_cluster_dist + 1)
            logger.info(
                f"Unfiltered insert clusters in contig {rname}: {str(len(clusters))}")
            # Get clusters
            jcss = self.getJuncClustersSeqs(
                clusters, rname, minseqs=mininsbreaks)
            logger.info(
                f"Insert clusters with cov > {str(mininsbreaks)} in contig {rname}: {str(len(jcss))}")
            # Extract reads and assemble consensus
            for jcs in jcss:
                # Check that extracted insert + flanking sequences are all
                # within len_threshold of the median length
                extr, coords = self.extractFlankingFromJcs(
                    jcs, rname, margin, 1.0-len_threshold, 1.0+len_threshold)
                aln = None
                if extr and coords:
                    # logging.debug(f"contig {rname} coordinates {str(coords[0])} {str(coords[1])}")
                    aln = self.spoaConsensusToFlanking(
                        extr, rname, coords[0], coords[1], margin=100, mode=1)
                if aln:
                    consseq, adjpos = findLongestInsert(
                        aln, rname, coords[0]-margin)
                    # Complain if predicted insert location from findLongestInsert
                    # is more than 5 bp away from the approximate insert location
                    problem = None
                    if adjpos < coords[0] and abs(adjpos - coords[0]) > 5:
                        logger.info(
                            f"Predicted insert further than 5 bp from preliminary position for {rname}, {str(adjpos)}")
                        problem = 1
                    elif adjpos > coords[1] and abs(adjpos - coords[1]) > 5:
                        logger.info(
                            f"Predicted insert further than 5 bp from preliminary position for {rname}, {str(adjpos)}")
                        problem = 1
                    if len(consseq) < 1:
                        logger.info(f"Consensus sequence of zero length for {rname}, {str(adjpos)}")
                        problem = 1
                    # Put together GFF entry
                    gffpos = adjpos + 1  # GFF convention
                    breakpointid = f"BREAK_POINTS_SUBREADS_{rname}_{str(gffpos)}_{str(len(consseq))}"
                    gfftype = "internal_eliminated_sequence_junction"
                    extrzmwcov = subreadNamesToZmwCoverage(
                        [r.id for r in extr])
                    attr = [f"ID={breakpointid}",
                            f"IES_length={str(len(consseq))}",
                            f"IES_subread_coverage={str(len(extr))}",
                            f"IES_zmw_coverage={str(extrzmwcov)}"]
                    # Get average coverage of region of interest
                    if self._alnformat == "bam":
                        regreads = self._alnfile.fetch(
                            str(rname), coords[0], coords[1]+1)
                        # plus one because coordinates are end-exclusive
                        regnames = [r.query_name for r in regreads
                                    if (not r.is_supplementary)
                                    and (not r.is_secondary)]
                        readcov = len(regnames)
                        zmwcov = subreadNamesToZmwCoverage(regnames)
                        attr.append(f"average_subread_coverage={str(readcov)}")
                        attr.append(f"average_zmw_coverage={str(zmwcov)}")
                    # Provisional approximate IES retention score
                    # R = IES+ / (IES+ + IES-)
                    # here the denominator is simply average coverage
                    provscore = None
                    if readcov:
                        if gfftype == "internal_eliminated_sequence_junction":
                            provscore = round(extrzmwcov/zmwcov, 4)
                    else:
                        logger.warn(
                            f"Region {rname} {str(coords[0])} {str(coords[1])} lacks read coverage")

                    consseq = SeqRecord(
                        Seq(consseq), id=breakpointid, description=";".join(attr)+";")
                    # Find pointers if present
                    if len(consseq) > 0:
                        gffpos, gffpos, pointerdict = self.reportAdjustPointers(
                            rname, gffpos, gffpos, consseq, breakpointid)
                        for i in pointerdict:
                            attr.append(f"{i}={str(pointerdict[i])}")

                    # Put together GFF entry
                    outarr = [str(rname),
                              "MILRAA",
                              gfftype,
                              str(gffpos),
                              str(gffpos),
                              str(provscore),
                              ".",
                              ".",
                              ";".join(attr)+";"]

                    if problem:
                        outseq_problem[consseq.id] = consseq
                        gff_problem.addEntry(outarr, None)
                    else:
                        outseq[consseq.id] = consseq
                        gff.addEntry(outarr, None)

                    counter += 1
                    if counter % 1000 == 0:
                        logger.info(f"Processed {counter} entries...")

        return(gff, outseq, gff_problem, outseq_problem)
