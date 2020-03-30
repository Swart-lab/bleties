#!/usr/bin/env python3

import re
import json
from collections import defaultdict

def getClips(cigar, pos):
    """Parse cigar string and alignment position and report left- and right-clipping
       Return list of tuples (start pos, end pos, and clip type)"""
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
    """Parse cigar string and alignment position and report insertions or deletions
       Return list of tuples (start pos, end pos, insert length, insertion or deletion)"""
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

    def __init__(self):
        """Constructor creates IesRecords, internaly represented by two dicts"""
        # dict to store counts of detected inserts keyed by evidence type
        # keys: contig -> startpos -> endpos -> ins_length -> evidence type -> count
        self._insDict = defaultdict(lambda: defaultdict( lambda: defaultdict ( lambda: defaultdict( lambda: defaultdict(int)))))

    def __str__(self):
        """String representation of IesRecords - as a JSON dump"""
        outstr = json.dumps(self._insDict, sort_keys = True, indent = 2)
        return(outstr)
