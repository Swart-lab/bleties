#!/usr/bin/env python3

import re
import json
from collections import defaultdict
import pysam

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

from bleties.SharedValues import SharedValues
from bleties.SharedFunctions import Gff

def getClips(cigar, pos):
    """Parse cigar string and alignment position and report left- and right-
    clipping. Return list of tuples (start pos, end pos, and clip type)

    Arguments:
    cigar -- Cigar string of alignment to process (str)
    pos -- Alignment start position on reference, 1-based numbering (int)
    """
    # Look for inserts that are on left or right ends of read (i.e. H or S clipping operations)
    outarr = []
    leftclipmatch = re.match(r"(\d+)[HS](\d+)M", cigar)
    if leftclipmatch:
        outarr.append((pos, pos, "HSM"))
    rightclipmatch = re.search(r"(\d+)M(\d+)[HS]$", cigar)
    if rightclipmatch:
        refconsumematches = re.findall(r"(\d+)[MDN\=X]", cigar) # find all ref-consuming operations
        rightclip_refconsumed = sum([int(refconsumematch) for refconsumematch in refconsumematches])
        outarr.append((pos + rightclip_refconsumed, pos+rightclip_refconsumed, "MHS"))
    return(outarr)

def getIndels(cigar, pos, minlength, qseq):
    """Parse cigar string and alignment position and report insertions or
    deletions. Return list of tuples (start pos, end pos, insert length,
    insertion or deletion, insert sequence). If it is a deletion then there is
    no insert sequence reported because that can be parsed from the reference.

    Arguments:
    cigar -- Cigar string of alignment to process (str)
    pos -- Alignment start position on reference, 1-based numbering (int)
    minlength -- Minimum length of indel to report (int)
    qseq -- Query sequence (str)
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
    Returns list of list (breakpoint, left junction seq, right junction seq)

    Arguments:
    iesgff -- Gff object listing insertions and deletions (output from reportPutativeIes)
    iesconsseq -- Dict of SeqRecords for indels (output from reportPutativeIes)
    refgenome -- Reference genome (dict of SeqRecord objects)
    flanklen -- Length of flanking sequence at junctions to report
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
        # If insertion, add 1 to end, because GFF convention is to record zero-
        # length features with start=end, and junction site is to right of
        # coordinate.
        if start == end: # insertion junction
            # end += 1
            indel = "I"
        elif start < end:
            indel = "D"
        else:
            raise Exception("Start cannot be less than end, feature " + breakpointid)
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
            flankrightseq = ref[ctg][end: end + flanklen].seq.lower() # TODO check if gff end inclusive
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
        # Put everything together and return
        outseqs.append([breakpointid,
                        str(flankleftseq),
                        str(flankrightseq),
                        str(iesconsseq[breakpointid].seq),
                        str(refunedited)
                        ])
    return(outseqs)


# TODO check flanks for potential pointers
def getPointers(iesgff,iesconsseq,ref):
    """Find potential pointer sequences at putative IES junctions

    Parameters
    ----------
    iesgff : Gff
        Gff object listing insertions and deletions, output from reportPutativeIes
    iesconsseq : dict
        dict of SeqRecords for indels (output from reportPutativeIes)
    ref : dict
        dict of SeqRecords, reference genome
    """
    pointdict = {}
    for breakpointid in iesgff:
        ctg = iesgff.getValue(breakpointid,'seqid')
        start = int(iesgff.getValue(breakpointid,'start'))
        end = int(iesgff.getValue(breakpointid,'end'))
        if start == end:
            indel = "I"
            iesseq = iesconsseq[breakpointid]
        elif start < end:
            indel = "D"
            iesseq = ref[ctg][start-1 : end]
        else:
            raise Exception(f"Start cannot be less than end, feature {breakpointid}")
        pointdict[breakpointid] = ""
        leftcheck = ""
        rightcheck = ""
        # check left of IES
        i = 0
        while iesseq[i] == ref[ctg][end+i]:
            leftcheck += iesseq[i]
            i += 1
            # break out of loop if the next position runs off edge
            if i >= len(iesseq) or end+i >= len(ref[ctg]):
                break
        # check right of IES
        i = -1
        if indel == "I":
            refstart = start
        elif indel == "D":
            refstart = start - 1
        while iesseq[i] == ref[ctg][refstart+i]:
            rightcheck = iesseq[i] + rightcheck
            i -= 1
            # break out of loop if next position runs off edge
            if -i > len(iesseq) or refstart+i < 0:
                break
        # take the longer putative pointer
        if len(leftcheck) < 2 and len(rightcheck) < 2:
            pointdict[breakpointid] = "none"
        elif len(leftcheck) > len(rightcheck) and len(leftcheck) >= 2:
            pointdict[breakpointid] = leftcheck
        elif len(rightcheck) > len(leftcheck) and len(rightcheck) >=2:
            pointdict[breakpointid] = rightcheck
        elif len(rightcheck) == len(leftcheck):
            pointdict[breakpointid] = "tie"
        else:
            logging.warn(f"Unexpected result in pointer search for breakpoint {breakpointid}")
    return(pointdict)

class IesRecords(object):
    """Records of putative IESs from mappings"""

    def __init__(self, alnfile, alnformat, refgenome):
        """Constructor creates IesRecords, internally represented by:
        _insDict -- dict to store counts of detected inserts/deletions, keyed
            by evidence type. Keys: contig (str) -> start pos (int) -> end
            pos (int) -> insert length (int) -> evidence type (str) -> count (int)
        _insSeqDict -- dict of sequences of detected inserts/deletions. Keys:
            contig (str) -> startpos (int) -> endpos (int) -> indel len (int) ->
            sequences (list of str)
        _alnfile -- as below
        _alnformat -- as below
        _refgenome -- as below

        Arguments:
        alnfile -- Alignment to parse (pysam.AlignmentFile)
        alnformat -- Format of the alignment, either "bam" or "sam"
        refgenome -- Reference genome sequences (Bio.SeqIO.SeqRecord)
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
        # pysam.AlignmentFile object representing the SAM/BAM mapping
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
        outstr = json.dumps(self._insDict, sort_keys = True, indent = 2)
        outstr_seq = json.dumps(self._insSeqDict, sort_keys=True, indent=2)
        return(outstr + "\n" + outstr_seq)

    def _addClipsFromCigar(self, rname, cigar, pos):
        """Check if alignment is clipped at the ends, and record the
        corresponding breakpoints in the _insDict dict, keyed by contig -> start
        position -> end position -> insert length -> evidence type, where
        evidence type is either "HSM" (left clip) or "MHS" (right clip).

        Arguments:
        rname -- Name of reference contig (str)
        cigar -- Cigar string of the current alignment record (str)
        pos -- Reference position of the current alignment record (int)
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

        Arguments:
        rname -- Name of reference contig (str)
        cigar -- Cigar string of the current alignment record (str)
        pos -- Reference position of the current alignment record (int)
        minlength -- Minimum length of indel for it to be recorded (int)
        qseq -- Query sequence of the read (str)
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
        Record them in the _insDict dict.

        Arguments:
        minlength -- Record only putative IESs of this length and above (int)
        """
        # parse CIGAR string
        for line in self._alnfile:
            if not line.flag & 4: # Skip unmapped reads
                if not line.flag & 256: # Skip secondary mappings
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

    def reportPutativeIes(self, mininsbreaks, mindelbreaks):
        """After clips and indels have been recorded, report putative IESs above
        the minimum coverage, and if the input alignment is a BAM file, then
        also report average coverage in the breakpoint region. Output is written
        to open file handle. Min coverage for deletion breakpoints is expected
        to be lower because we are mapping to somatic genome in the typical use
        case, and reads with alternative excisions are thought to be rare.

        Return:
        SharedFunctions.Gff object
        Dict of consensus sequences for putative IESs

        Arguments:
        mininsbreaks -- Minimum breakpoint coverage to report potential insertion (int)
        mindelbreaks -- Minimum breakpoint coverage to report potential deletion (int)
        gff -- SharedFunctions.Gff object
        """
        # Create lists to hold SeqRecord objects and GFF output
        gff = Gff()
        outseq = {}
        # Parse the dict and report putative IESs above min coverage
        # We only check breakpoints which are completely spanned by a read ("I" or "D" operations)
        # however we also report supporting counts from HSM and MSH type mappings
        # TODO Reduce code duplication here, incorporate Gff module
        for ctg in sorted(self._insDict):
            for ins_start in sorted(self._insDict[ctg]):
                for ins_end in sorted(self._insDict[ctg][ins_start]):
                    for ins_len in sorted(self._insDict[ctg][ins_start][ins_end]):
                        for evidencetype in self._insDict[ctg][ins_start][ins_end][ins_len]:
                            countvalue = self._insDict[ctg][ins_start][ins_end][ins_len][evidencetype]
                            # If the breakpoint is an insert type
                            if evidencetype == "I":
                                if countvalue >= mininsbreaks:
                                    breakpointid = "_".join(["BREAK_POINTS",str(ctg),str(ins_start),str(ins_end),str(ins_len)])
                                    # Prepare attributes list of key-value pairs
                                    attr = ["ID="+breakpointid,
                                            "IES_length="+str(ins_len)
                                           ]
                                    attr.append("cigar=I "+str(countvalue))
                                    # Extend the cigar evidence counts with clipped reads that are also at this junction
                                    # Clippings are recorded with nominal insert length of zero
                                    attr.extend(["cigar="+cigartype+" "+str(self._insDict[ctg][ins_start][ins_end][0][cigartype])
                                                 for cigartype in sorted(self._insDict[ctg][ins_start][ins_end][0])
                                                ])
                                    # Get read coverage from BAM file; SAM does not allow random access
                                    if self._alnformat == "bam":
                                        readcov = self._alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_end)) # TODO: Check for off-by-one errors
                                        attr.append("average_coverage="+str(readcov))
                                    # Get indel consensus
                                    consseq = self.reportIndelConsensusSeq(ctg, ins_start, ins_end, ins_len)
                                    consseq.id = breakpointid
                                    consseq.description = ";".join(attr)+";"
                                    outseq[consseq.id] = consseq
                                    # Build GFF entry
                                    outarr = [str(ctg),            # 1 seqid
                                              "MILRAA",            # 2 source
                                              "junction",          # 3 type
                                              str(ins_start),      # 4 start
                                              str(ins_end),        # 5 end
                                              str(countvalue),     # 6 score - in this case, breakpoint counts for insert operation only
                                              ".",                 # 7 strand
                                              ".",                 # 8 phase
                                              ";".join(attr)+";"   # 9 attributes
                                              ]
                                    gff.addEntry(outarr, None)
                            # If the breakpoint is a deletion type
                            elif evidencetype == "D":
                                if countvalue >= mindelbreaks:
                                    # Add 1 because both start and end are inclusive
                                    del_len = int(ins_end) - int(ins_start) + 1 # TODO: Check for off-by-one errors
                                    breakpointid = "_".join(["BREAK_POINTS",str(ctg),str(ins_start),str(ins_end),str(del_len)])
                                    # Build attributes field
                                    attr = ["ID="+breakpointid,
                                            "IES_length="+str(del_len)
                                           ]
                                    attr.append("cigar=D "+str(countvalue))
                                    # Look for left- and rightclips that fall on the deletion boundaries
                                    if self._insDict[ctg][ins_start][ins_start][ins_len].get("MHS") and int(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]) > 0:
                                        attr.append("cigar=MHS "+str(self._insDict[ctg][ins_start][ins_start][ins_len]["MHS"]))
                                    if self._insDict[ctg][ins_end][ins_end][ins_len].get("HSM") and int(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]) > 0:
                                        attr.append("cigar=HSM "+str(self._insDict[ctg][ins_end][ins_end][ins_len]["HSM"]))
                                    # Add coverage value
                                    if self._alnformat == "bam":
                                        readcov = self._alnfile.count(str(ctg), start=int(ins_start)-1, stop=int(ins_end))
                                        attr.append("average_coverage="+str(readcov))
                                    # Get indel consensus
                                    consseq = self.reportIndelConsensusSeq(ctg, ins_start, ins_end, del_len)
                                    consseq.id = breakpointid
                                    consseq.description = ";".join(attr)+";"
                                    outseq[consseq.id] = consseq
                                    # Build GFF entry
                                    outarr = [str(ctg),            # 1 seqid
                                              "MILRAA",            # 2 source
                                              "region",            # 3 type
                                              str(ins_start),      # 4 start
                                              str(ins_end),        # 5 end
                                              str(countvalue),     # 6 score - in this case, breakpoint counts for delete operation only
                                              ".",                 # 7 strand
                                              ".",                 # 8 phase
                                              ";".join(attr)+";"   # 9 attributes
                                              ]
                                    gff.addEntry(outarr, None)
        return(gff, outseq)

    def reportIndelConsensusSeq(self, ctg, indelstart, indelend, indellen):
        """Report consensus of indel sequence
        Returns: Bio.SeqRecord object

        Arguments:
        ctg -- name of contig (str)
        indelstart -- start position, 1-based (int)
        indelend -- end position, 1-based (int)
        indellen -- length of indel (int)
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

    def reportIndelReadMismatchPc(self, ctg, indelstart, indelend, indellen):
        """Report sequence mismatch % of query reads containing indel at a
        specific position vs. reads without indel.
        This is to flag indels that may originate from paralogs and hence are
        probably not true IESs.

        Arguments:
        ctg -- name of contig (str)
        indelstart -- start position of indel, 1-based (int)
        indelend -- end position, 1 based (int)
        indellen -- length of indel (int)

        Returns:
        ins_mm -- list of floats, mismatch % of query reads containing indel at
                  target position
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
