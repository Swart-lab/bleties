#!/usr/bin/env python3

import logging

from bleties.SharedFunctions import *
from bleties.SharedValues import SharedValues


logger = logging.getLogger("Milroc")


class IesCorrelationsByRead(object):
    def __init__(self, gfffile, alnfile):
        """Construct object to store IES correlation data.
        
        Parameters
        ----------
        gfffile : str
            Path to GFF3 file of IES annotations
        alnfile : pysam.AlignmentFile
            Alignment file containing mapping of interest
        """
        self._gfffile = gfffile
        self._alnfile = alnfile
        self._gff = Gff()
        self._gff.file2gff(gfffile)
        # Rekey the GFF entries by contig and coordinates
        self._gffByContig = defaultdict(
                lambda: defaultdict(
                    lambda: defaultdict(
                        list)))
        for gffid in self._gff:
            ctg = self._gff.getValue(gffid, "seqid")
            start = int(self._gff.getValue(gffid, "start"))
            end = int(self._gff.getValue(gffid, "end"))
            self._gffByContig[ctg][start][end].append(gffid)
        # Initialize dict to record IESs observed per read
        self._perRead = defaultdict( # read ID
                lambda: defaultdict(
                    list)) # list of IES IDs
        # Initialize dict to record IES co-occurrences
        self._perIes = defaultdict( # IES ID
                lambda: defaultdict( # ID of other IESs that co-occur with first
                    int)) # Number of co-occurrences/co-observations


    def countIesCooccurrences(self):
        """Iterate through alignment and count IES co-occurrences per read
        """
        # Iterate across all reads
        for alnrec in self._alnfile.fetch():
            ctg = alnrec.reference_name
            start = int(alnrec.reference_start) + 1 # convert to 1-based inclusive
            end = int(alnrec.reference_end) # 1-based inclusive
            qname = alnrec.query_name
            if not qname:
                logger.warn("Query in alignment without a query_name field")
            # Record info on each aligned read
            self._perRead[qname]['ref'] = ctg
            self._perRead[qname]['start'] = start
            self._perRead[qname]['end'] = end
            # Check if this contig has IESs defined
            if ctg in self._gffByContig: 
                for iesstart in self._gffByContig[ctg]:
                    if iesstart >= start and iesstart <= end:
                        for iesend in self._gffByContig[ctg][iesstart]:
                            if iesend <= end:
                                res = getOperationAtRefPos(iesstart,
                                        start,
                                        alnrec.cigarstring,
                                        1, 1)
                                if res:
                                    op, op_len = res
                                    if op == "I":
                                        # Read has insert corresponding to IES position
                                        for iesid in self._gffByContig[ctg][iesstart][iesend]:
                                            self._perRead[qname]['present'].append(iesid)
                                        # TODO account for insert length
                                        # TODO also work for deletions
                                    elif op == "M":
                                        # No IES in the defined position
                                        for iesid in self._gffByContig[ctg][iesstart][iesend]:
                                            self._perRead[qname]['absent'].append(iesid)
            # For all the IESs present on this read, tally co-occurrences with
            # other IESs observed
            for ies1 in self._perRead[qname]['present']:
                for ies2 in self._perRead[qname]['present']:
                    self._perIes[ies1][ies2] += 1


    def dump(self, filename):
        """Dump internal data to JSON for troubleshooting

        Parameters
        ----------
        filename : str
            Path to file to write JSON
        """
        import json
        with open(filename, "w") as fh:
            fh.write(json.dumps({"_perRead" : self._perRead,
                "_perIes" : self._perIes}, indent=4))


    def summarizePerRead(self):
        """Summarize statistics per read for plotting

        Returns
        -------
        list
            list of lists containing summary IES statistics per read.
            Fields: read name, reference contig, start, end, number of IESs
            with inserts present, number of IESs with inserts absent
        """
        out = []
        for qname in self._perRead:
            out.append([qname,
             self._perRead[qname]['ref'],
             self._perRead[qname]['start'],
             self._perRead[qname]['end'],
             len(self._perRead[qname]['present']),
             len(self._perRead[qname]['absent'])])
        return(out)

    # TODO: ? Graphical representation of IES co-occurrence
    # TODO: Flag IESs that have co-occurrences more or less than expected
