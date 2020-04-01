#!/usr/bin/env python3

import re
import json
from collections import defaultdict
import pysam

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

def getIndels(cigar, pos, minlength):
    """Parse cigar string and alignment position and report insertions or
    deletions. Return list of tuples (start pos, end pos, insert length,
    insertion or deletion)

    Arguments:
    cigar -- Cigar string of alignment to process (str)
    pos -- Alignment start position on reference, 1-based numbering (int)
    minlength -- Minimum length of indel to report (int)
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
                outarr.append((ins_pos_start, ins_pos_start, ins_len, "I"))
        # We also look for delete operations above a min length
        # These are already present in the reference, so the "insert length" is 0
        if cigmatch.group(2) == "D": # If delete operation,
            del_len = int(cigmatch.group(1)) # Length of current deletion
            if del_len >= minlength:
                # Get start and end pos of deletion
                del_pos_start = pos + ref_consumed # 1-based, insert starts to right of position
                del_pos_end = pos + ref_consumed + del_len # 1-based
                outarr.append((del_pos_start, del_pos_end, 0, "D"))
        # Count ref and query consumed _after_ the insert has been accounted for
        if cigmatch.group(2) in REFCONSUMING:
            ref_consumed += int(cigmatch.group(1))
        if cigmatch.group(2) in QUERYCONSUMING:
            que_consumed += int(cigmatch.group(1))
    return(outarr)

class IesRecords(object):
    """Records of putative IESs from mappings"""

    def __init__(self, alnfile, alnformat):
        """Constructor creates IesRecords, internally represented by a dict

        Arguments:
        alnfile -- Alignment to parse (pysam.AlignmentFile)
        alnformat -- Format of the alignment, either "bam" or "sam"
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
        self._alnfile = alnfile
        self._alnformat = alnformat

    def __str__(self):
        """Report summary stats of IesRecords object"""
        dictlen = len(self._insDict)
        nref = self._alnfile.nreferences
        mapped = self._alnfile.mapped
        return("bleties.IesRecords object with dict of length "
                + str(dictlen)
                + " and alignment with " 
                + str(nref) 
                + " references and "
                + str(mapped) 
                + " mapped reads")

    def dump(self):
        """Data dump of IesRecords._insDict in JSON format"""
        outstr = json.dumps(self._insDict, sort_keys = True, indent = 2)
        return(outstr)

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

    def _addIndelsFromCigar(self, rname, cigar, pos, minlength):
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
        """
        # Look for inserts that are completely spanned by the read (i.e. I operations)
        indelarr = getIndels(cigar, pos, minlength)
        if len(indelarr) > 0:
            for (indelstart, indelend, indellen, indeltype) in indelarr:
                self._insDict[rname][indelstart][indelend][indellen][indeltype] += 1

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

            # Find left and right clips and record them
            self._addClipsFromCigar(rname, line.cigarstring, pos)
            # Find indels (putative IESs) over the minimum length and record them
            self._addIndelsFromCigar(rname, line.cigarstring, pos, minlength)
            # # if int(total_mismatch) - int(total_i) < 0: 
            #     # Sanity check - mismatches include inserts, but cannot be fewer than inserts
            #     # print ("Uh-oh!")

    def reportPutativeIes(self, minbreaks, fh):
        """After clips and indels have been recorded, report putative IESs above
        the minimum coverage, and if the input alignment is a BAM file, then
        also report average coverage in the breakpoint region. Output is written
        to open file handle.

        Arguments:
        minbreaks -- Minimum breakpoint coverage to report (int)
        fh -- File handle for writing output
        """
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
                                if countvalue >= minbreaks:
                                    # Prepare attributes list of key-value pairs
                                    attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start)+"_"+str(ins_end),
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
                                    fh.write("\t".join(outarr)+"\n")
                            # If the breakpoint is a deletion type
                            elif evidencetype == "D":
                                if countvalue >= minbreaks:
                                    del_len = int(ins_end) - int(ins_start) # TODO: Check for off-by-one errors
                                    attr = ["ID=BREAK_POINTS_"+str(ctg)+"_"+str(ins_start)+"_"+str(ins_end),
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
                                    fh.write("\t".join(outarr)+"\n")
