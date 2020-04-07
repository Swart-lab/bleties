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
    # Lists of reference- and query-consuming operations
    REFCONSUMING = ['M', 'D', 'N', '=', 'X']
    QUERYCONSUMING = ['M', 'I', 'S', '=', 'X']
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
                ins_pos_start = pos + ref_consumed # 1-based, insert is to the right of position
                # Get the sequence of the insert
                ins_seq = qseq[que_consumed:que_consumed + ins_len] # 0-based, following pysam convention
                outarr.append((ins_pos_start, ins_pos_start, ins_len, "I", ins_seq))
        # We also look for delete operations above a min length
        # These are already present in the reference, so the "insert length" is 0
        if cigmatch.group(2) == "D": # If delete operation,
            del_len = int(cigmatch.group(1)) # Length of current deletion
            if del_len >= minlength:
                # Get start and end pos of deletion
                del_pos_start = pos + ref_consumed # 1-based, insert starts to right of position
                del_pos_end = pos + ref_consumed + del_len # 1-based
                outarr.append((del_pos_start, del_pos_end, 0, "D", "")) # If deletion, no insert sequence reported
        # Count ref and query consumed _after_ the insert has been accounted for
        if cigmatch.group(2) in REFCONSUMING:
            ref_consumed += int(cigmatch.group(1))
        if cigmatch.group(2) in QUERYCONSUMING:
            que_consumed += int(cigmatch.group(1))
    return(outarr)

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
                    if self._insDict[rname][indelstart][indelend][0]["D"]:
                        # Record the sequence to the insSelfDict
                        # Convert from 1-based to 0-based numbering for slicing
                        indellen = int(indelend) - int(indelstart)
                        indelseq = str(refctgseq[int(indelstart)-1:int(indelend)-1]) # TODO: Check for off-by-one errors
                        self._insSeqDict[rname][indelstart][indelend][indellen].append(indelseq)

    def findPutativeIes(self, minlength):
        """Search alignment for clips and indels to identify putative IESs.
        Record them in the _insDict dict. 

        Arguments:
        minlength -- Record only putative IESs of this length and above (int)
        """
        # parse CIGAR string
        for line in self._alnfile:
            pos = int(line.reference_start) + 1 # Convert from 0-based numbering in pysam to 1-based in GFF3 and SAM
            rname = line.reference_name # Get reference name
            total_mismatch = line.get_tag("NM") # Get number of mismatches
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

        Arguments:
        mininsbreaks -- Minimum breakpoint coverage to report potential insertion (int)
        mindelbreaks -- Minimum breakpoint coverage to report potential deletion (int)
        """
        # Create lists to hold SeqRecord objects and GFF output
        outseq = []
        outgff = []
        # Parse the dict and report putative IESs above min coverage
        # We only check breakpoints which are completely spanned by a read ("I" or "D" operations)
        # however we also report supporting counts from HSM and MSH type mappings
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
                                    outseq.append(consseq)
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
                                    outgff.append(outarr)
                            # If the breakpoint is a deletion type
                            elif evidencetype == "D":
                                if countvalue >= mindelbreaks:
                                    del_len = int(ins_end) - int(ins_start) # TODO: Check for off-by-one errors
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
                                    outseq.append(consseq)
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
                                    outgff.append(outarr)
        return(outgff, outseq)

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

