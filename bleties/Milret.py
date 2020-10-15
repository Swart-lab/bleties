#!/usr/bin/env python3

import re
from collections import defaultdict
import logging

from bleties.SharedValues import SharedValues
from bleties.SharedFunctions import *


# Define logger
logger = logging.getLogger("Milret")


class IesRetentionsMacOnly(object):
    def __init__(self, gfffile, alnfile):
        """Construct object to store IES retention data.
        Using only mapping to MAC genome assembly. Unlike ParTIES MIRET which
        uses both mapping to somatic and germline genomes

        Parameters
        ----------
        gfffile : str
            Path to GFF3 file of IES annotations (str)
        alnfile : pysam.AlignmentFile
            Alignment file containing mapping of interest
        """
        self._gfffile = gfffile
        self._alnfile = alnfile
        # Read in GFF file to memory
        self._gff = Gff()
        self._gff.file2gff(gfffile)
        # Initialize dict to hold counts of operations at IES junctions
        self._countsDict = defaultdict(  # ID of GFF feature
            lambda: defaultdict(    # Operation (M, D, I)
                list)               # list of op lengths at junction
        )
        # Initialize dict to hold retention scores calculated from counts
        self._scoresDict = defaultdict(float)

    def dump(self, dumpfile):
        """Dump internal objects to JSON file for troubleshooting

        Parameters
        ----------
        dumpfile : str
            Path to output file to write JSON
        """
        import json
        with open(dumpfile, "w") as fh:
            fh.write(json.dumps({"countsDict": self._countsDict,
                                 "scoresDict": self._scoresDict},
                                indent=4))

    def findMappingOps(self):
        """Find mapping operations at the IES junctions, and count how many of
        each type.

        Match operations (M) that span the IES junction are treated as
        representing IES- form, because the mapping reference is somatic
        genome. Insert operations (I) that are exactly at the IES junction are
        treated as representing IES+ form.
        """
        # TODO: Set cutoff for minimum length of match oepration, and/or min
        # distance for the match boundaries from the IES junction
        for gffid in self._gff._gffDict:  # Each IES ID
            seqid = self._gff.getValue(gffid, 'seqid')
            start = int(self._gff.getValue(gffid, 'start'))
            end = int(self._gff.getValue(gffid, 'end'))
            # GFF allows zero-length features. However, to fetched alignments
            # from pysam.AlignmentFile, we need nonzero interval
            if start == end:
                end += 1
            # Get alignments that span this position
            # fetch() method uses 1-based SAM coordinates, unlike rest of Pysam!
            for alnrec in self._alnfile.fetch(seqid, start, end):
                # For this aligned read, which alignment operation spans this
                # position?
                # Convert 0-based pysam numbering to 1-based in GFF 
                # TODO check off-by-one errors
                res = getOperationAtRefPos(start,
                                           alnrec.reference_start+1,
                                           alnrec.cigarstring,
                                           1,  # TODO: Let user choose thresholds
                                           1)
                if res:  # If there is no operation, will return None
                    op, op_len = res
                    # Record the count
                    self._countsDict[gffid][op].append(op_len)

    def calculateRetentionScores(self):
        """Calculate retention scores from counts of I and M operations per site
        after findMappingOps() has been applied.

        Equation: R = IES+ / (IES+ + IES-)
        """

        for gffid in self._gff:
            # Only calculate for junction features, i.e. inserts
            if self._gff.getValue(gffid, "start") == self._gff.getValue(gffid, "end"):
                if gffid in self._countsDict:
                    iesplus = 0
                    iesminus = 0
                    if 'M' in self._countsDict[gffid]:
                        iesminus = len(self._countsDict[gffid]['M'])
                    if 'I' in self._countsDict[gffid]:
                        iesplus = len(self._countsDict[gffid]['I'])
                    if iesplus + iesminus > 0:
                        score = round(float(iesplus / (iesplus + iesminus)), 4)
                    else:
                        score = None
                    self._scoresDict[gffid] = score
            else:
                logger.debug(
                    f"Ignored GFF entry {gffid} when calculating retention scores, not a junction")

    def calculateRetentionScoresMatchLengths(self, threshold=0.05):
        """Calculate retention scores from counts of I and M operations per site
        after findMappingOps() has been applied.

        Only count inserts that match the length of reported IES in the input
        GFF file, within a defined threshold. This assumes that the input GFF
        file is produced by MILRAA and contains an `IES_length` field in the
        attributes.

        Equation: R = IES+ / (IES+ + IES-)

        In the output table, the retention score will account for IES length,
        but the total number of "I" operations will be unchanged (includes both
        inserts that match the length and those that do not).
        """

        for gffid in self._gff:
            # Only calculate for junction features, i.e. inserts
            if self._gff.getValue(gffid, "start") == self._gff.getValue(gffid, "end"):
                # Get defined IES length
                ieslength = self._gff.getAttr(gffid, "IES_length")
                if ieslength:  # None is returned if attribute absent
                    ieslength = mean_of_number_list(ieslength)
                if gffid in self._countsDict:
                    iesplus = 0
                    iesminus = 0
                    if 'M' in self._countsDict[gffid]:
                        iesminus = len(self._countsDict[gffid]['M'])
                    if 'I' in self._countsDict[gffid]:
                        iesplus_lengthmatched = [i for i in self._countsDict[gffid]['I']
                                                 if i >= ieslength*(1-threshold)
                                                 and i <= ieslength*(1+threshold)]
                        iesplus = len(iesplus_lengthmatched)
                    if iesplus + iesminus > 0:
                        score = round(float(iesplus / (iesplus + iesminus)), 4)
                    else:
                        score = None
                    self._scoresDict[gffid] = score
            else:
                logger.debug(
                    f"Ignored GFF entry {gffid} when calculating retention scores, not a junction")

    def reportRetentionScores(self, outfile):
        """Report retention scores after running calculateRetentionScores().
        Writes to filehandle.

        Parameters
        ----------
        outfile : str
            Path to file to write results
        """
        # Create header line
        headerarr = ['ID', 'score']
        headerarr.extend(SharedValues.ALLCIGAROPS)
        with open(outfile, "w") as fh:
            fh.write("\t".join(headerarr)+"\n")
            # Report for each IES junction
            for gffid in sorted(self._scoresDict):
                # Report scores
                outarr = [gffid, str(self._scoresDict[gffid])]
                # Report counts per CIGAR op, zero if not recorded for this junction
                for op in SharedValues.ALLCIGAROPS:
                    if self._countsDict[gffid][op]:
                        outarr.append(str(len(self._countsDict[gffid][op])))
                    else:
                        outarr.append("0")
                fh.write("\t".join(outarr) + "\n")
