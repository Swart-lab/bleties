#!/usr/bin/env python3

import re
import json
from collections import defaultdict
import pysam
import logging

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


def getClips(cigar, pos):
    """Parse cigar string and alignment position and report left- and right-
    clipping.

    Parameters
    ----------
    cigar : str
        CIGAR string of alignment to process
    pos : int
        Alignment start position on reference, 1-based numbering

    Returns
    -------
    list
        list of tuples (int, int, str) representing (start position, end
        position, clipping type)
    """
    # Look for inserts that are on left or right ends of read (i.e. H or S clipping operations)
    outarr = []
    leftclipmatch = re.match(r"(\d+)[HS](\d+)M", cigar)
    if leftclipmatch:
        outarr.append((pos, pos, "HSM"))
    rightclipmatch = re.search(r"(\d+)M(\d+)[HS]$", cigar)
    if rightclipmatch:
        refconsumematches = re.findall(r"(\d+)[MDN\=X]", cigar) # find all ref-consuming operations
        rightclip_refconsumed = sum([int(refconsumematch) 
            for refconsumematch in refconsumematches])
        outarr.append((pos + rightclip_refconsumed, 
            pos+rightclip_refconsumed, 
            "MHS"))
    return(outarr)


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
    outarr = [] # Array to store tuples of results
    # Initialize counters for how much query or reference seq is consumed
    ref_consumed = 0
    que_consumed = 0
    # Split cigar string into individual operations
    cigs = re.findall(r"\d+[\w\=]", cigar)
    for cig in cigs:
        cigmatch = re.match(r"(\d+)([\w\=])",cig) # Get number and operation
        # We want to look for insert operations above a min length,
        # map their positions on the reference, and count how many reads
        # support a given putative insert
        # This requires that we count the operations that consume reference
        # and add it to the POS field
        if cigmatch.group(2) == "I": # If insert operation,
            ins_len = int(cigmatch.group(1)) # Length of the current insert
            if ins_len >= minlength: # Check that insert is above min length
                # Get the start and end positions of the insert
                ins_pos_start = pos + ref_consumed - 1 # 1-based, insert is to the right of position
                # Get the sequence of the insert
                ins_seq = qseq[que_consumed:que_consumed + ins_len] # 0-based, following pysam convention
                outarr.append((ins_pos_start, ins_pos_start, ins_len, "I", ins_seq))
        # We also look for delete operations above a min length
        # These are already present in the reference, so the "insert length" is 0
        if cigmatch.group(2) == "D": # If delete operation,
            del_len = int(cigmatch.group(1)) # Length of current deletion
            if del_len >= minlength:
                # Get start and end pos of deletion
                del_pos_start = pos + ref_consumed # 1-based, deletion starts ON this position
                del_pos_end = pos + ref_consumed + del_len - 1 # 1-based, deletion ends ON this position
                outarr.append((del_pos_start, del_pos_end, 0, "D", "")) # If deletion, no insert sequence reported
        # Count ref and query consumed _after_ the insert has been accounted for
        if cigmatch.group(2) in SharedValues.REFCONSUMING:
            ref_consumed += int(cigmatch.group(1))
        if cigmatch.group(2) in SharedValues.QUERYCONSUMING:
            que_consumed += int(cigmatch.group(1))
    return(outarr)


def getIndelJunctionSeqs(iesgff,iesconsseq,ref,flanklen):
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
        ctg = iesgff.getValue(breakpointid,'seqid')
        start = int(iesgff.getValue(breakpointid,'start'))
        end = int(iesgff.getValue(breakpointid,'end'))
        flankleftseq = ""
        flankrightseq = ""
        indel = ""
        refunedited = ""
        # Check if this is insertion junction or deletion region
        # GFF convention is to record zero-length features with start=end, and
        # junction site is to right of coordinate.
        if start == end: # insertion junction
            indel = "I"
        elif start < end:
            indel = "D"
        else:
            raise Exception(f"Start cannot be less than end, feature {breakpointid}")
        # Get the left flanking junction on reference
        # Check whether sequence with flanking will run off the end
        if start >= flanklen: # TODO Check off-by-one
            if indel == "I": # For insert, feature starts on RIGHT of coordinate
                flankleftseq = ref[ctg][start - flanklen: start].seq.lower()
            elif indel == "D": # For deletion, feature starts ON the coordinate
                flankleftseq = ref[ctg][start-1-flanklen:start-1].seq.lower()
        else: # Pad the left side
            flankleftseq = (flanklen - start) * "-" + ref[ctg][0:start].seq.lower()
        # Get the right flanking junction on reference
        # Check whether sequence with flanking will run off the end
        if (end + flanklen) <= len(ref[ctg]):
            flankrightseq = ref[ctg][end: end + flanklen].seq.lower()
        else:
            flankrightseq = ref[ctg][end:].seq.lower()
        # Add the IES sequence fragment
        if indel == "I":
            flankleftseq = flankleftseq + iesconsseq[breakpointid][0:flanklen].seq.upper()
            flankrightseq = iesconsseq[breakpointid][-flanklen:].seq.upper() + flankrightseq
        elif indel == "D":
            # In GFF, start and end are 1-based and BOTH inclusive
            flankleftseq = flankleftseq + ref[ctg][start-1: start-1+flanklen].seq.upper()
            flankrightseq = ref[ctg][end-flanklen:end].seq.upper() + flankrightseq
        # Get the reference sequence
        if indel =="I":
            refunedited = ref[ctg][start-flanklen:end+flanklen].seq.lower()
        elif indel =="D":
            # Start minus one because for deletion, start is ON the deleted region
            refunedited = ref[ctg][start-flanklen-1:end+flanklen].seq.lower()
        # Find pointers if present
        (pointer, pointerstart, pointerend)  = getPointers(ref[ctg], start, end, iesconsseq[breakpointid], breakpointid)
        if pointerstart != start:
            logger.info(f"Position of pointer {breakpointid} has been adjusted")
            start = pointerstart
            end = pointerend
        # Convert pointers to TA junctions if possible
        (tastart, taend, tapointer) = adjustPointerTA(pointerstart, pointerend, pointer)
        # Maximize pointer lengths
        (ppstart, ppend, pppointer) = adjustPointerMaxlength(ref[ctg], start, end, pointer, iesconsseq[breakpointid])
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
        Sequence of contig with putative IES junction/region
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
        iesseq = seq[start-1 : end]
    else:
        raise Exception(f"Start cannot be less than end, feature {name}")

    # Remove gap characters from ies sequence
    if isinstance(iesseq, str):
        iesseq = iesseq.replace("-","")
    elif isinstance(iesseq, SeqRecord):
        iesseq = iesseq.seq.ungap("-")
    else:
        raise Exception(f"iesseq must be of type str or SeqRecord but is {type(iesseq)}")

    pointer = ""
    pointerstart, pointerend = start, end
    leftcheck = ""
    rightcheck = ""
    # check left of IES
    i = 0
    while iesseq[i] == seq[end+i]:
        # print(f"checkleft {name} {seq} {str(i)} {iesseq[i]} {seq[end+i]}") # for debugging
        # print(iesseq)
        leftcheck += iesseq[i]
        i += 1
        # break out of loop if the next position runs off edge
        if i >= len(iesseq) or end+i >= len(seq):
            break
    # check right of IES
    i = -1
    if indel == "I": # beacuse GFF puts zero-length features to right of coordinate
        refstart = start
    elif indel == "D":
        refstart = start - 1
    while iesseq[i] == seq[refstart+i]:
        # print(f"checkright {name} {seq} {str(i)} {iesseq[i]} {seq[refstart+i]}") # for debugging
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
    elif len(rightcheck) > len(leftcheck) and len(rightcheck) >=2:
        pointer = rightcheck
        # If the insert position is to the right of the putative pointer seq,
        # report adjusted pointer position such that the insert position always
        # lies to the left of the putative pointer seq on the MDS.
        pointerstart = pointerstart - len(pointer)
        pointerend = pointerend - len(pointer)
    elif len(rightcheck) == len(leftcheck):
        logger.info(f"Breakpoint {name} has potential pointers on both sides")
        pointer = "tie" # if both sides could potentially have a pointer
    else:
        logger.warn(f"Unexpected result in pointer search for breakpoint {name}")
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
        iesseq = seq[start - 1 : end]
    elif start == end:
        indel = "I"
        if not iesseq:
            raise Exception ("IES sequence is missing")
    # Look upstream of the pointer to try and extend match
    i = -1
    if indel == "I": # beacuse GFF puts zero-length features to right of coordinate
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


def alignSeqsMuscle(seqlist, muscle_path="muscle"):
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
    elif isinstance(seqlist[0], Bio.Seq.Seq):
        seqrecs = [SeqRecord(i) for i in seqlist]
    elif isinstance(seqlist[0], Bio.SeqRecord.SeqRecord):
        seqrecs = seqlist
    else:
        raise Exception("sequence list must comprise str, Seq, or SeqRecord objects")
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
        raise Exception("sequence list must comprise str, Seq, or SeqRecord objects")
    aln = MultipleSeqAlignment(seqrecs)
    alninf = AlignInfo.SummaryInfo(aln)
    alncons = alninf.dumb_consensus(threshold=threshold)
    alnconsrec = SeqRecord(alncons)
    return(alnconsrec)
    


class IesRecords(object):
    """Records of putative IESs from mappings"""


    def __init__(self, alnfile, alnformat, refgenome):
        """Constructor for IesRecords

        Internally represented by:
        _insDict -- dict to store counts of detected inserts/deletions, keyed
            by evidence type. Keys: contig (str) -> start pos (int) -> end
            pos (int) -> insert length (int) -> evidence type (str) -> count (int)
        _insSeqDict -- dict of sequences of detected inserts/deletions. Keys:
            contig (str) -> startpos (int) -> endpos (int) -> indel len (int) ->
            sequences (list of str)
        _alnfile -- as below
        _alnformat -- as below
        _refgenome -- as below

        Parameters
        ----------
        alnfile : pysam.AlignmentFile
            Alignment to parse
        alnformat : str
            Format of the alignment, either "bam" or "sam"
        refgenome : Bio.SeqIO.SeqRecord
            Reference genome sequences
        """
        # dict to store counts of detected inserts keyed by evidence type
        # keys: contig -> startpos -> endpos -> insert length -> evidence type -> count
        self._insDict = defaultdict(              # contig
                lambda: defaultdict(              # startpos
                    lambda: defaultdict (         # endpos
                        lambda: defaultdict(      # insert length
                            lambda: defaultdict(  # evidence type
                                int)              # count
                            )
                        )
                    )
                )
        # dict to store sequences of detected inserts/deletions
        self._insSeqDict = defaultdict(      # contig
                lambda: defaultdict(         # startpos
                    lambda: defaultdict(     # endpos
                        lambda: defaultdict( # indel length
                            list)            # list of sequences (str)
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
        insdictlen = len(self._insDict)
        insseqdictlen = len(self._insSeqDict)
        nref = self._alnfile.nreferences
        alnformat = self._alnformat
        mapped = self._alnfile.mapped
        return("bleties.IesRecords object with "
                + " insDict of length "
                + str(insdictlen)
                + " and insSeqDict of length "
                + str(insseqdictlen)
                + " and alignment of format "
                + str(alnformat)
                + " with "
                + str(nref)
                + " references and "
                + str(mapped)
                + " mapped reads")


    def dump(self):
        """Data dump of IesRecords._insDict in JSON format"""
        outstr = json.dumps({"_insDict" : self._insDict, "_insSeqDict" : self._insSeqDict}, sort_keys = True, indent = 2)
        return(outstr)


    def _addClipsFromCigar(self, rname, cigar, pos):
        """Check if alignment is clipped at the ends, and record the
        corresponding breakpoints in the _insDict dict, keyed by contig -> start
        position -> end position -> insert length -> evidence type, where
        evidence type is either "HSM" (left clip) or "MHS" (right clip).

        Parameters
        ----------
        rname : str
            Name of reference contig
        cigar : str
            CIGAR string of the current alignment record
        pos : int
            Reference position of the current alignment record
        """
        # Look for inserts that are on left or right read ends (i.e. H or S clipping operations)
        cliparr = getClips(cigar, pos)
        if len(cliparr) > 0:
            for (clipstart, clipend, cliptype) in cliparr:
                self._insDict[rname][clipstart][clipend][0][cliptype] += 1


    def _addIndelsFromCigar(self, rname, cigar, pos, minlength, qseq):
        """Check if alignment contains indels above minimum length, and record
        the corresponding breakpoints relative to the reference, and the insert
        length. If the indel is an insert, insert length > 0. If the indel is a
        deletion, insert length = 0. Recorded in the _insDict dict, keyed by
        contig -> start pos -> end pos -> insert length -> evidence type, where
        evidence type is either "I" (insertion) or "D" (deletion).

        Parameters
        ----------
        rname : str
            Name of reference contig
        cigar : str
            Cigar string of the current alignment record
        pos : int
            Reference position of the current alignment record
        minlength : int
            Minimum length of indel for it to be recorded
        qseq : str
            Query sequence of the read
        """
        # Look for inserts that are completely spanned by the read (i.e. I operations)
        indelarr = getIndels(cigar, pos, minlength, qseq)
        if len(indelarr) > 0:
            for (indelstart, indelend, indellen, indeltype, indelseq) in indelarr:
                self._insDict[rname][indelstart][indelend][indellen][indeltype] += 1
                # If insertion, record the inserted sequence to dict
                if indeltype == "I":
                    self._insSeqDict[rname][indelstart][indelend][indellen].append(indelseq)


    def _getDeletedSequences(self):
        """Record sequences of deletions. Sequences of insertions are recorded
        when parsing each alignment, because they are extracted from the query.
        However deletions are extracted from the reference, so they can be done
        in a single pass.

        Returns: Changes object in-place
        """
        for rname in self._insDict: # Get contig name
            if rname not in self._refgenome:
                logger.error(f"Sequence {rname} not found in reference genome")
            # Get reference sequence from Fasta file
            refctgseq = self._refgenome[rname].seq
            # For each start and stop position
            for indelstart in self._insDict[rname]:
                for indelend in self._insDict[rname][indelstart]:
                    # If there is a deletion
                    if 0 in self._insDict[rname][indelstart][indelend]:
                        if "D" in self._insDict[rname][indelstart][indelend][0]:
                            # Record the sequence to the insSelfDict
                            # Convert from 1-based to 0-based numbering for slicing
                            # Plus one because end and start in GFF features are BOTH inclusive
                            indellen = int(indelend) - int(indelstart) + 1
                            # Get sequence of indel sequence, note that end is also inclusive
                            indelseq = str(refctgseq[int(indelstart)-1:int(indelend)]) # TODO: Check for off-by-one errors
                            self._insSeqDict[rname][indelstart][indelend][indellen].append(indelseq)


    def findPutativeIes(self, minlength):
        """Search alignment for clips and indels to identify putative IESs.
        Record them in the _insDict dict and sequences in _insSeqDict

        Parameters
        ----------
        minlength : int
            Record only putative IESs of this length and above
        """
        # parse CIGAR string
        for line in self._alnfile:
            if (not line.is_unmapped) and (not line.is_secondary) and (not line.is_supplementary):
                pos = int(line.reference_start) + 1 # Convert from 0-based numbering in pysam to 1-based in GFF3 and SAM
                rname = line.reference_name # Get reference name
                # total_mismatch = line.get_tag("NM") # Get number of mismatches
                # total_i = 0
                qseq = line.query_sequence # Get query sequence
                # Find left and right clips and record them
                self._addClipsFromCigar(rname, line.cigarstring, pos)
                # Find indels (putative IESs) over the minimum length and record them
                self._addIndelsFromCigar(rname, line.cigarstring, pos, minlength, qseq)
                # # if int(total_mismatch) - int(total_i) < 0:
                #     # Sanity check - mismatches include inserts, but cannot be fewer than inserts
                #     # print ("Uh-oh!")

        # Go through insDict and extract sequences of deleted regions, too
        self._getDeletedSequences()


    def reportPutativeIesInsertFuzzy(self, mininsbreaks, mindelbreaks, dist_threshold=0.05):
        """Report putative IESs where insert lengths do not match exactly

        Parameters
        ----------
        mininsbreaks
        mindelbreaks
            Same as reportPutativeIes()
        dist_threshold : float
            Distance threshold for clustering, in range (0,1)
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
                    ins_seqs.extend(dd[i])
                clust_seqs, clust_ids = get_clusters_from_seqlist(ins_seqs, dist_threshold)
                # For each cluster, report a putative IES
                for clust in clust_seqs:
                    counts = defaultdict(int)
                    for i in clust:
                        counts[len(i)] += 1
                    totalcount = sum(counts.values())
                    maxcounts = [l for l in counts if counts[l] == max(counts.values())]

                    # Report only if above minimum
                    if totalcount >= mininsbreaks:
                        # Put together ID for this breakpoint
                        if len(counts.keys()) == 1: # cluster of one length
                            prefix = f"BREAK_POINTS"
                        elif len(counts.keys()) > 1:
                            prefix = f"BREAK_POINTS_FUZZY"
                        breakpointid = "_".join([prefix, str(ctg), str(ins_start), str(ins_end)] + [str(l) for l in counts.keys()])
                        gfftype = "internal_eliminated_sequence_junction"

                        # Attributes list of key-value pairs
                        # report modal IES length
                        attr = ["ID=" + breakpointid,
                                "IES_length="+"_".join([str(l) for l in maxcounts])]
                        # Report number of counts per insert length
                        attr.append("cigar=" + " ".join([str(l) + "I*" + str(counts[l]) for l in counts]))

                        # Get indel consensus
                        if len(counts.keys()) == 1: # cluster of one length
                            # take dumb consensus since all seqs are same length
                            consseq = alnDumbFromSeqs(clust)
                        else:
                            # Otherwise use Muscle (potentially slower) to align
                            consseq = alnFromSeqs(clust)

            elif ins_end > ins_start:
                # Deletion feature - need to get count from _insDict
                # This is otherwise same as reportPutativeIes
                del_len = ins_end - ins_start + 1 # GFF is end-inclusive
                totalcount = self._insDict[ctg][ins_start][ins_end][0]["D"]
                if totalcount >= mindelbreaks:
                    breakpointid = "_".join(["BREAK_POINTS",str(ctg),str(ins_start),str(ins_end),str(del_len)])
                    gfftype = "internal_eliminated_sequence"
                    consseq = SeqRecord(Seq(dd[del_len][0], generic_dna))
                    # Build attributes field
                    attr = ["ID="+breakpointid,
                            "IES_length="+str(del_len)]
                    attr.append("cigar=" + str(del_len) + "D*" + str(totalcount))
                    # Look for left- and rightclips that fall on the deletion boundaries
                    # if self._insDict[ctg][ins_start][ins_start][ins_len].get("MHS") and int(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]) > 0:
                    #     attr.append("cigar=MHS*" + str(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]))
                    # if self._insDict[ctg][ins_end][ins_end][ins_len].get("HSM") and int(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]) > 0:
                    #     attr.append("cigar=HSM*" + str(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]))

            else:
                raise Exception("feature cannot start after it ends")

            # Record to Gff
            if breakpointid:
                # Get average coverage of region of interest
                if self._alnformat == "bam":
                    readcov = self._alnfile.count(str(ctg),
                            start=int(ins_start) - 1,
                            stop=int(ins_end))
                    attr.append("average_coverage="+str(readcov))

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
                (pointer, pointerstart, pointerend)  = getPointers(self._refgenome[ctg], ins_start, ins_end, consseq, breakpointid)
                if pointerstart != ins_start:
                    logger.info(f"Position of pointer {breakpointid} has been adjusted")
                    ins_start = pointerstart
                    ins_end = pointerend
                if pointer: # Add pointer seq to attributes field if present
                    attr.append("pointer_seq=" + str(pointer))
                # Convert pointers to TA junctions if possible
                (tastart, taend, tapointer) = adjustPointerTA(pointerstart, pointerend, pointer)
                # Add TA pointer sequence to attributes field if present
                if tapointer:
                    attr.extend(["ta_pointer_seq=" + str(tapointer),
                        "ta_pointer_start=" + str(tastart),
                        "ta_pointer_end=" + str(taend)])
                # Maximize pointer lengths
                (ppstart, ppend, pppointer) = adjustPointerMaxlength(self._refgenome[ctg], ins_start, ins_end, pointer, consseq)
                if pppointer:
                    attr.extend(["pp_pointer_seq=" + str(pppointer),
                        "pp_pointer_start=" + str(ppstart),
                        "pp_pointer_end=" + str(ppend)])

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
        # Inserts: fuzzy cluster start/end pos (start == end), then cluster
        # ins_lengths.
        # Deletions: fuzzy cluster start and end positions separately.
        # This is somewhat trickier, because of the way the data are structured


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
        # however we also report supporting counts from HSM and MSH type mappings
        # TODO Reduce code duplication here
        for rec in nested_dict_to_list(self._insDict):
            [ctg, ins_start, ins_end, ins_len, evidencetype, countvalue] = rec
            if ctg not in self._refgenome:
                logger.error(f"Sequence {ctg} not found in reference genome")

            breakpointid = None
            attr = []
            indel_len = 0

            # If the breakpoint is an insert type
            if evidencetype == "I" and countvalue >= mininsbreaks:
                breakpointid = "_".join(["BREAK_POINTS",str(ctg),str(ins_start),str(ins_end),str(ins_len)])
                indel_len = ins_len
                gfftype = "internal_eliminated_sequence_junction"
                # Prepare attributes list of key-value pairs
                attr = ["ID="+breakpointid,
                        "IES_length="+str(ins_len)]
                attr.append("cigar=" + str(ins_len) + "I*" + str(countvalue))
                # Extend the cigar evidence counts with clipped reads that are also at this junction
                # Clippings are recorded with nominal insert length of zero
                attr.extend(["cigar=" + cigartype + "*"+str(self._insDict[ctg][ins_start][ins_end][0][cigartype])
                             for cigartype in sorted(self._insDict[ctg][ins_start][ins_end][0])
                            ])

            # If the breakpoint is a deletion type
            elif evidencetype == "D" and countvalue >= mindelbreaks:
                # Add 1 because both start and end are inclusive
                del_len = int(ins_end) - int(ins_start) + 1 # TODO: Check for off-by-one errors
                indel_len = del_len
                breakpointid = "_".join(["BREAK_POINTS",str(ctg),str(ins_start),str(ins_end),str(del_len)])
                gfftype = "internal_eliminated_sequence"
                # Build attributes field
                attr = ["ID="+breakpointid,
                        "IES_length="+str(del_len)]
                attr.append("cigar=" + str(del_len) + "D*" + str(countvalue))
                # Look for left- and rightclips that fall on the deletion boundaries
                if self._insDict[ctg][ins_start][ins_start][ins_len].get("MHS") and int(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]) > 0:
                    attr.append("cigar=MHS*" + str(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]))
                if self._insDict[ctg][ins_end][ins_end][ins_len].get("HSM") and int(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]) > 0:
                    attr.append("cigar=HSM*" + str(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]))

            # If the breakpoint has been defined, report it to the GFF file
            if breakpointid and attr:
                # Get read coverage from BAM file; SAM does not allow random access
                if self._alnformat == "bam":
                    readcov = self._alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_end)) # TODO: Check for off-by-one errors
                    attr.append("average_coverage="+str(readcov))

                # Provisional approximate IES retention score
                # R = IES+ / (IES+ + IES-)
                # here the denominator is simply average coverage
                if readcov:
                    if gfftype == "internal_eliminated_sequence_junction":
                        provscore = round(countvalue/readcov, 4)
                    elif gfftype == "internal_eliminated_sequence":
                        provscore = round((readcov-countvalue)/readcov, 4)

                # Get indel consensus
                consseq = self.reportIndelConsensusSeq(ctg, ins_start, ins_end, indel_len)
                consseq.id = breakpointid
                consseq.description = ";".join(attr)+";"
                outseq[consseq.id] = consseq

                # Find pointers if present
                (pointer, pointerstart, pointerend)  = getPointers(self._refgenome[ctg], ins_start, ins_end, consseq, breakpointid)
                if pointerstart != ins_start:
                    logger.info(f"Position of pointer {breakpointid} has been adjusted")
                    ins_start = pointerstart
                    ins_end = pointerend
                if pointer: # Add pointer seq to attributes field if present
                    attr.append("pointer_seq=" + str(pointer))
                # Convert pointers to TA junctions if possible
                (tastart, taend, tapointer) = adjustPointerTA(pointerstart, pointerend, pointer)
                # Add TA pointer sequence to attributes field if present
                if tapointer:
                    attr.extend(["ta_pointer_seq=" + str(tapointer),
                        "ta_pointer_start=" + str(tastart),
                        "ta_pointer_end=" + str(taend)])
                # Maximize pointer lengths
                (ppstart, ppend, pppointer) = adjustPointerMaxlength(self._refgenome[ctg], ins_start, ins_end, pointer, consseq)
                if pppointer:
                    attr.extend(["pp_pointer_seq=" + str(pppointer),
                        "pp_pointer_start=" + str(ppstart),
                        "pp_pointer_end=" + str(ppend)])

                # Build GFF entry
                outarr = [str(ctg),            # 1 seqid
                          "MILRAA",            # 2 source
                          gfftype,             # 3 type
                          str(ins_start),      # 4 start
                          str(ins_end),        # 5 end
                          str(provscore),      # 6 score - provisional IES retention score
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

        # From list of sequences as str, make a list of SeqRecord objects
        seqrecs = [SeqRecord(Seq(i, generic_dna)) for i in self._insSeqDict[ctg][indelstart][indelend][indellen]]
        # Make a pseudo-alignment from the list of SeqRecord objects
        aln = MultipleSeqAlignment(seqrecs)
        # Summarize alignment, and get dumb consensus sequence
        alninf = AlignInfo.SummaryInfo(aln)
        alncons = alninf.dumb_consensus()
        # alnconsrec = SeqRecord(alncons, id=consname)
        alnconsrec = SeqRecord(alncons)
        return(alnconsrec)


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
                    [SeqRecord(Seq(i, generic_dna))
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



    def reportIndelReadMismatchPc(self, ctg, indelstart, indelend, indellen):
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

        Returns
        -------
        list, list
            Two lists:
            ins_mm -- list of floats, mismatch % of query reads containing indel
                    at target position
            non_mm -- list of floats, mismatch % of query reads without indel at pos
        """
        MIN_IES_LEN = 10                                                        # TODO: replace magic number
        # GFF allows start==end, but pysam does not recognise
        dummyend = indelend
        if indelend == indelstart:
            dummyend += 1
        # Get segments that overlap indel of interest
        itrr = self._alnfile.fetch(ctg, indelstart - 1, dummyend - 1) # minus 1 for pysam uses 0-based coords
        segs = [seg for seg in itrr] # Segments aligning to position of interest
        # initialize lists to report mismatch percentages
        non_mm = []
        ins_mm = []
        # for each segment, check if it contains indel at position of interest
        for seg in segs:
            # if seg.query_sequence is None:
            #     print(f"Query {seg.query_name} is None!")
            indels = getIndels(seg.cigarstring, int(seg.reference_start), MIN_IES_LEN, seg.query_sequence)
            indelcoords = set([(int(indel[0]),int(indel[1])) for indel in indels])
            mismatch_pc = 100 * float(seg.get_tag("NM")) / float(seg.query_length) # number of mismatchs / query length * 100 pc
            if (indelstart - 1, indelend - 1) in indelcoords:
                ins_mm.append(mismatch_pc)
            else:
                non_mm.append(mismatch_pc)
        return(ins_mm, non_mm)


