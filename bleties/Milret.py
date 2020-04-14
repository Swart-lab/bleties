#!/usr/bin/env python3

import re
from collections import defaultdict
from bleties.SharedValues import SharedValues
from bleties.SharedFunctions import SharedFunctions
from bleties.SharedFunctions import Gff

def getOperationAtRefPos(reftargetpos, refstartpos, cigar, mininslength, minmatchlength):
    """Identify whether a given reference position is covered by a reference-
    consuming operation. 
    Returns the operation that covers that position. If no operation covers
    that position, nothing is returned

    Arguments:
    reftargetpos -- Target position on the reference sequence, 1-based (int)
    refstartpos -- Start position of the alignment on reference, 1-based (int)
    cigar -- CIGAR string of the alignment (str)
    mininslength -- Minimum length of insert to report (int)
    minmatchlength - Minimum length of match to report (int)
    """
    
    curr_int_start = refstartpos
    curr_int_end = refstartpos
    # Split cigar string into individual operations
    cigs = re.findall(r"\d+[\w\=]", cigar)
    for cig in cigs:
        cigmatch = re.match(r"(\d+)([\w\=])",cig) # Get number and operation
        # If reference is consumed:
        if cigmatch.group(2) in SharedValues.REFCONSUMING:
            # Update the current interval
            curr_int_start = curr_int_end
            curr_int_end = curr_int_end + int(cigmatch.group(1))
        # Otherwise if query is consumed
        else:
            curr_int_start = curr_int_end # current operation has extent zero on ref

        # If current operation is ref-consuming
        if curr_int_end > curr_int_start:
            # Check whether the target position is contained in the current interval
            if reftargetpos in range(curr_int_start,curr_int_end): # TODO check off-by-one errors
                if int(cigmatch.group(1)) > minmatchlength:
                    return(cigmatch.group(2))
        # If current operation interval is zero (i.e. not ref-consuming)
        elif curr_int_end == curr_int_start:
            # and it matches exactly the target poosition
            if reftargetpos == curr_int_end:
                if int(cigmatch.group(1)) > mininslength:
                    return(cigmatch.group(2))

class IesRetentionsMacOnly(object):
    def __init__(self, gfffile, alnfile):
        """Construct object to store IES retention data.
        Using only mapping to MAC genome assembly. Unlike ParTIES MIRET which
        uses both mapping to somatic and germline genomes

        Arguments:
        gfffile - Path to GFF3 file of IES annotations (str)
        alnfile - pysam.AlignmentFile object containing mapping of interest
        """
        self._gfffile = gfffile
        self._alnfile = alnfile
        # Read in GFF file to memory
        self._gff = Gff()
        self._gff.file2gff(gfffile)
        # Initialize dict to hold counts of operations at IES junctions
        self._countsDict = defaultdict( # ID of GFF feature
                lambda: defaultdict(    # Operation (M, D, I)
                    int)                # Count of op at junction
                )
        # Initialize dict to hold retention scores calculated from counts
        self._scoresDict = defaultdict(float) 

    def findMappingOps(self):
        """Find mapping operations at the IES junctions, and count how many of
        each type. Match operations (M) that span the IES junction are treated
        as representing IES- form, because the mapping reference is somatic
        genome. Insert operations (I) that are exactly at the IES junction are
        treated as representing IES+ form.
        """
        # TODO: Set cutoff for minimum length of match oepration, and/or min
        # distance for the match boundaries from the IES junction
        # TODO: Only count insert operations that are within X bases in length
        # different from the reported IES length
        for gffid in self._gff._gffDict: # Each IES ID
            seqid = self._gff._gffDict[gffid]['seqid']
            start = int(self._gff._gffDict[gffid]['start'])
            end = int(self._gff._gffDict[gffid]['end'])
            # GFF allows zero-length features. However, to fetched alignments
            # from pysam.AlignmentFile, we need nonzero interval
            if start == end:
                end += 1
            # TODO: Also get average coverage at this position, using pysam
            # Get alignments that span this position
            for alnrec in self._alnfile.fetch(seqid, start, end):
                # For this aligned read, which alignment operation spans this
                # position?
                res = getOperationAtRefPos(start,
                                           alnrec.reference_start+1, # Convert 0-based pysam numbering to 1-based in GFF TODO check off-by-one errors
                                           alnrec.cigarstring,
                                           1, # TODO: Let user choose thresholds
                                           1)
                if res: # If there is no operation, will return None
                    # Record the count
                    self._countsDict[gffid][res] +=1 
    
    def calculateRetentionScores(self):
        """Calculate retention scores from counts of I and M operations per site
        after findMappingOps() has been applied.
        Equation: R = IES+ / (IES+ + IES-)
        """

        for gffid in self._countsDict:
            iesplus = 0
            iesminus = 0
            if 'M' in self._countsDict[gffid]:
                iesminus  = self._countsDict[gffid]['M']
            if 'I' in self._countsDict[gffid]:
                iesplus = self._countsDict[gffid]['I']
            if iesplus + iesminus > 0:
                score = iesplus / (iesplus + iesminus)
            else:
                score = None
            self._scoresDict[gffid] = score

    def reportRetentionScores(self, fh):
        """Report retention scores after running calculateRetentionScores().
        Writes to filehandle.

        Arguments:
        fh - Filehandle to write results
        """
        # Create header line
        headerarr = ['ID', 'score']
        headerarr.extend(SharedValues.ALLCIGAROPS)
        fh.write("\t".join(headerarr)+"\n")
        # Report for each IES junction
        for gffid in sorted(self._scoresDict):
            # Report scores
            outarr = [gffid, str(self._scoresDict[gffid])]
            # Report counts per CIGAR op, zero if not recorded for this junction
            for op in SharedValues.ALLCIGAROPS:
                if self._countsDict[gffid][op]:
                    outarr.append(str(self._countsDict[gffid][op]))
                else:
                    outarr.append("0")
            fh.write("\t".join(outarr) + "\n")
