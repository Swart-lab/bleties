#!/usr/bin/env python3

"""Functions shared between all Bleties scripts"""

import re
from collections import defaultdict
from operator import itemgetter
from bleties.SharedValues import SharedValues

class SharedFunctions():
    def returntrue():
        return (True)

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

    def addEntry(self, linearr, gffid):
        """ Add single GFF entry to Gff object

        Arguments:
        linearr -- List of GFF3 fields, length 9 (list)
        gffid -- ID attribute of the entry
        """
        # Check that GFF line has only nine fields
        if len(linearr) != 9:
            raise Exception('GFF3 input encountered with incorrect number of fields')
        idval = ""
        idsearch = re.search(r"ID=([^;]+);", linearr[8])
        if idsearch:
            idval = idsearch.group(1)
        elif gffid:
            idval = gffid
        else:
            raise Exception('GFF3 input without ID field')
        # Make a dict of each GFF3 column name and value
        self._gffDict[idval] = dict(zip(SharedValues.GFF3COLUMNS, linearr))
        # Split attributes field into key-value pairs and make a dict of it
        attrs = re.findall(r"([^;=]+)=([^;=]+)[;$]", linearr[8])
        attrkeys = [match[0] for match in attrs]
        attrvals = [match[1] for match in attrs]
        attrdict = dict(zip(attrkeys, attrvals))
        self._gffDict[idval]['attrdict'] = attrdict

    def getValue(self, gffid, column):
        """Get value of column for a given GFF entry.
        GFF entry is specified by the ID.
        Returns the value.

        Arguments:
        gffid -- ID of the GFF entry
        column -- Name of the column to retrieve
        """
        if column in SharedValues.GFF3COLUMNS:
            if gffid in self._gffDict:
                return(self._gffDict[gffid][column])
            else:
                raise Exception("Unknown GFF3 ID " + gffid)
        else:
            raise Exception("Unknown GFF3 column name " + column)

    def getAttr(self, gffid, attribute):
        """Get value from attributes field of a GFF entry.

        Arguments:
        gffid -- ID of the GFF entry
        attribute -- Key of the attribute requested
        """
        if gffid in self._gffDict:
            if attribute in self._gffDict[gffid]['attrdict']:
                return(self._gffDict[gffid]['attrdict'][attribute])
            else:
                raise Exception("Unknown attribute " + attribute + "for GFF3 ID " + gffid)
        else:
            raise Exception("Unknown GFF3 ID " + gffid)


    def list2gff(self, gfflist):
        """Add entries to Gff object from a list of GFF3 lines.
        
        Arguments:
        gfflist -- List of GFF3 records (list of str)
        """
        # Iterate through list, split each line into columns
        for line in gfflist:
            if not re.match(r"#", line): # Skip header lines
                line = line.rstrip() # Strip trailing whitespace
                linearr = line.split("\t")
                self.addEntry(linearr, None)

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

    def file2gff(self, filename):
        """Read in GFF3 file directly to Gff object
        """
        with open(filename, "r") as fh:
            slurp = [line.rstrip() for line in fh]
            self.list2gff(slurp)

    def gff2file(self, filename):
        """Write Gff object directly to GFF3 file
        """
        fh = open(filename,"w")
        outarr = self.gff2list()
        for line in outarr:
            fh.write(line+"\n")
        fh.close()

    def gff2fh(self, fh):
        """Write Gff object to open filehandle
        """
        outarr = self.gff2list()
        for line in outarr:
            fh.write(line+"\n")
