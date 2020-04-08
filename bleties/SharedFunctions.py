#!/usr/bin/env python3

"""Functions shared between all Bleties scripts"""

import re
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


