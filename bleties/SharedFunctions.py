#!/usr/bin/env python3

"""Functions shared between all Bleties scripts"""

import re
from collections import defaultdict
from operator import itemgetter
from bleties.SharedValues import SharedValues

class SharedFunctions():
    def isMatchAtRefPos(reftargetpos, refstartpos, cigar, mininslength, minmatchlength):
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

class Gff(object):
    def __init__(self):
        """Construct Gff object
        Gff object is a dict keyed by ID field of each GFF entry.
        Each dict is keyed by the GFF column name, e.g. 'strand', 'start'.
        There is an additional key 'attrdict', which is a dict of key-value pairs
        in the 'attributes' field of the GFF column. This assumes that the 
        attributes keys are unique.
        """
        # Create dict 
        self._gffDict = defaultdict(dict)

    def list2gff(self, gfflist):
        """Add entries to Gff object from a list of GFF3 lines.
        
        Arguments:
        gfflist -- List of GFF3 records (list of str)
        """
        # Iterate through list, split each line into columns
        for line in gfflist:
            line = line.rstrip() # Strip trailing whitespace
            linearr = line.split("\t")
            # Check that GFF line has only nine fields
            if len(linearr) != 9:
                raise Exception('GFF3 input encountered with incorrect number of fields')
            # Get ID from attr field
            idsearch = re.search(r"ID=([^;]+);", linearr[8])
            idval = idsearch.group(1)
            # Make a dict of each GFF3 column name and value
            self._gffDict[idval] = dict(zip(SharedValues.GFF3COLUMNS, linearr))
            # Split attributes field into key-value pairs and make a dict of it
            attrs = re.findall(r"([^;=]+)=([^;=]+)[;$]", linearr[8])
            attrkeys = [match[0] for match in attrs]
            attrvals = [match[1] for match in attrs]
            attrdict = dict(zip(attrkeys, attrvals))
            self._gffDict[idval]['attrdict'] = attrdict

    def gff2list(self):
        """Write Gff3 object to list of strings, for printing.
        Returns list of str representing GFF3 entries, sorted.
        Effectively reconstitutes the input to list2gff()
        """

        outarr=[]
        for gffid in self._gffDict: # ID value
            linearr = [self._gffDict[gffid][colname] for colname in SharedValues.GFF3COLUMNS]
            # Make numeric fields int - necessary for correct sorting later
            linearr[3] = int(linearr[3])
            linearr[4] = int(linearr[4])
            outarr.append(linearr)
        # Sort by seqid, start, and end fields in that order
        outarr = sorted(outarr, key=itemgetter(0,3,4))
        # Covnert fields back to str and join back into a tab-separated string
        outarr = [map(str, linearr) for linearr in outarr]
        outarr = ["\t".join(linearr) for linearr in outarr]
        return(outarr)
