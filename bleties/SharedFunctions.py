#!/usr/bin/env python3

"""Functions shared between all Bleties scripts"""

import re
from collections import defaultdict
from operator import itemgetter
from Bio.Align import PairwiseAligner
from Bio.Cluster import treecluster
from numpy import array
import logging

from bleties.SharedValues import SharedValues


logger = logging.getLogger("SharedFunctions")


class SharedFunctions():
    def returntrue():
        return (True)


def nested_dict_to_list(d):
    """Flatten a nested dict of dicts to lists of key-key-value sets

    Parameters
    ----------
    d : dict
        dict of dicts (of dicts ... ) with a nested structure, e.g.
        key1 -> key2 -> ... -> val

    Returns
    -------
    list
        list where each element is a list of the key value sets, e.g.
        [key1, key2, ..., val]
    """
    out = []

    def recc(dd, p=[]):
        if not (type(dd) is dict or type(dd) is defaultdict):
            out.append(p + [dd])
            return(p + [dd])
        else:
            for k in dd:
                recc(dd[k], p + [k])
    recc(d)
    return(out)


def nested_dict_to_list_fixed_depth(d, depth):
    """Flatten a nested dict of dicts to lists of key-key-value sets

    Parameters
    ----------
    d : dict
        dict of dicts (of dicts ... ) with a nested structure, e.g.
        key1 -> key2 -> ... -> val
    depth : int
        Depth to which to flatten the nested dict

    Returns
    -------
    list
        list where each element is a list of the key value sets, e.g.
        [key1, key2, ..., dict], but only to the depth specified. E.g. depth = 1
        would yield [key1, dict]
    """
    out = []

    def recc(dd, p=[]):
        if len(p) == depth or not (type(dd) is dict or type(dd) is defaultdict):
            out.append(p + [dd])
            return(p + [dd])
        else:
            for k in dd:
                recc(dd[k], p + [k])
    recc(d)
    return(out)


def get_clusters(in_list, cluster_type: str, width: int):
    """Cluster a list of integers into groups not more than _width_ apart

    Note that this will also report clusters of length 1

    Parameters
    ----------
    in_list : list
        List of integers
    cluster_type : str
        Type of clustering to perform, either "bp" (by base pair width) or "pc"
        (by percentage difference)
    width : int
        Maximum distance (exclusive) between _adjacent_ members of a cluster,
        i.e. width=2 would mean that 101 and 102 would be clustered together.
        width=1 would mean no clustering.

    Returns:
    list
        list of lists, each sub-list containing the values of each cluster
    """
    # convert to ints just in case
    in_list = [int(x) for x in in_list]
    in_list = sorted(in_list)
    merged = []
    hold = []
    prev = None
    for i in in_list:
        if prev:
            # if distance from previous value less than condition
            # append to the list
            condition = False
            if cluster_type == "bp":
                condition = (i - prev < width)
            elif cluster_type == "pc":
                condition = within_percent(i, prev, width)
            if condition:
                hold.append(prev)
            else:
                # else make cluster and start new cluster
                hold.append(prev)
                merged.append(hold)
                hold = []
        prev = i
    if prev:  # sweet up last entry
        hold.append(prev)
        merged.append(hold)
    return(merged)


def within_percent(num1, num2, percent: int):
    """Compare two numeric values by percentage difference

    Return True if they are mutually within x-percent of each other

    Parameters
    ----------
    num1 : int or float
    num2 : int or float
    percent : int
        Percentage difference between the two. Mutual difference!

    Returns
    -------
    bool
        True if num1 and num2 are within percent of each other
    """
    # Sort numerically, convert to float just in case
    compsorted = sorted([float(num1), float(num2)])
    lower = 1 - (float(percent)/100)
    upper = 1 + (float(percent)/100)
    if compsorted[0] * upper > compsorted[1] * lower:
        return(True)
    else:
        return(False)


def get_clusters_from_seqlist(seqlist, dist_threshold=0.05):
    """Cluster a list of sequences by a distance identity threshold

    Parameters
    ----------
    seqlist : list
        list of sequences as str
    dist_threshold : float
        Max distance value to retain, branches above this length in the 
        hierarchical clustering tree will be cut.

    Returns
    -------
    list
        list of lists - input sequences now grouped by cluster
    list
        list of int - cluster memberships of the originally input list
    """
    if len(seqlist) == 1:
        # Skip alignment if there is only one sequence
        return([seqlist], [0])
    else:
        aligner = PairwiseAligner()
        aligner.mode = "local"

        # Convert sequence list to distance matrix
        distmatrix = []
        for seq1 in seqlist:
            row = []
            for seq2 in seqlist:
                maxlen = max([len(seq1), len(seq2)])
                # Take percentage identity of pairwise alignment score (match base
                # +1, all other operations +0) over the longer sequence in pair
                idval = aligner.align(seq1, seq2).score / maxlen
                distval = 1 - idval  # convert to distance fraction
                row.append(distval)
            distmatrix.append(row)
        # Hierarchical clustering from the distance matrix
        htree = treecluster(data=None, distancematrix=array(distmatrix))
        # Find number of branches with length longer than threshold, and add 1
        # to get number of cuts
        cuts = 1 + len([htree[i].distance for i in range(len(htree))
                        if htree[i].distance > dist_threshold])
        clust_ids = list(htree.cut(cuts))
        clust_seqs_dict = defaultdict(list)
        for i in range(len(seqlist)):
            clust_seqs_dict[clust_ids[i]] += [seqlist[i]]
        # Convert dict of lists to list of lists
        clust_seqs = [clust_seqs_dict[i] for i in clust_seqs_dict]
        return(clust_seqs, clust_ids)


def getOperationAtRefPos(reftargetpos, refstartpos, cigar, mininslength, minmatchlength):
    """Identify whether a given reference position is covered by a reference-
    consuming operation.

    Parameters
    ----------
    reftargetpos : int
        Target position on the reference sequence, 1-based
    refstartpos : int
        Start position of the alignment on reference, 1-based
    cigar : str
        CIGAR string of the alignment
    mininslength : int
        Minimum length of insert to report
    minmatchlength : int
        Minimum length of match to report

    Returns
    -------
    str
        Operation that covers that position. If no operation covers that
        position, nothing is returned
    int
        Length of the operation.
    """
    curr_int_start = refstartpos - \
        1  # the minus-one is necessary to get this to work, TODO figure out why!
    curr_int_end = refstartpos - 1
    # print(f"{str(refstartpos)} {str(reftargetpos)} {cigar}") # diagnostic mode
    # Split cigar string into individual operations
    cigs = re.findall(r"\d+[\w\=]", cigar)
    for cig in cigs:
        cigmatch = re.match(r"(\d+)([\w\=])", cig)  # Get number and operation
        # If reference is consumed:
        if cigmatch.group(2) in SharedValues.REFCONSUMING:
            # Update the current interval
            curr_int_start = curr_int_end
            curr_int_end = curr_int_end + int(cigmatch.group(1))
        # Otherwise if query is consumed
        else:
            curr_int_start = curr_int_end  # current operation has extent zero on ref

        # If current operation is ref-consuming
        if curr_int_end > curr_int_start:
            # Check whether the target position is contained in the current interval
            # TODO check off-by-one errors
            if reftargetpos in range(curr_int_start, curr_int_end):
                if int(cigmatch.group(1)) > minmatchlength:
                    # print(f"{str(curr_int_start)} {str(curr_int_end)} {str(reftargetpos)} {cig}") # diagnostic mode
                    return(cigmatch.group(2), int(cigmatch.group(1)))
        # If current operation interval is zero (i.e. not ref-consuming)
        elif curr_int_end == curr_int_start:
            # and it matches exactly the target poosition
            if reftargetpos == curr_int_end:
                if int(cigmatch.group(1)) > mininslength:
                    # print(f"{str(curr_int_start)} {str(curr_int_end)} {str(reftargetpos)} {cig}") # diagnostic mode
                    return(cigmatch.group(2), int(cigmatch.group(1)))


def getCigarOpQuerySeqs(qseq, cigartuples, rstart, target_op="S"):
    """Get sequence segments from a query sequence that correspond to a
    specific CIGAR operation

    Parameters
    ----------
    qseq : str
        Query sequence
    cigartuples : list
        List of tuples of ints (operation, operation length) where operation is
        CIGAR operation numbered 0-9 by BAM convention. Output from
        pysam.AlignedSegment.cigartuples
    rstart : int
        Reference start position, 0-based
    target_op : str
        Name of CIGAR operation for which to report the sequence(s)

    Returns
    -------
    list
        list of tuples (str, int, int, int, int) corresponding to: (sequence
        segment, query start, query end, ref start, ref end) of each segment
        corresponding to a specific CIGAR operation. Coordinates are 0-based.
    """
    op2bam = {"M": 0, "I": 1, "D": 2, "N": 3,
              "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    if target_op not in op2bam:
        raise Exception(f"Operation {target_op} not a valid CIGAR operation")
    QUERY_CONSUMING = [0, 1, 4, 7, 8]
    REF_CONSUMING = [0, 2, 3, 7, 8]

    # initialize placeholders
    rend = rstart
    qstart = 0  # 0-based coordinates
    qend = qstart

    out = []
    for op, oplen in cigartuples:
        if op in QUERY_CONSUMING:
            qend = qstart + oplen
        if op in REF_CONSUMING:
            rend = rstart + oplen

        if op == op2bam[target_op]:
            out.append((qseq[qstart:qend],
                        qstart, qend,
                        rstart, rend))

        # Update counters
        qstart = qend
        rstart = rend

    if len(out) > 0:
        return(out)


def mean_of_number_list(numbers, delim="_"):
    """Report arithmetic mean of a string of integers that are separated by a
    common delimiter. Used here to get mean of IES lengths reported in MILRAA
    output GFF attribute.

    Parameters
    ----------
    numbers : str
        List of numbers separated by a common delimiter
    delim : str
        Character used as delimiter

    Returns
    -------
    int
        Arithmetic mean rounded to nearest integer
    """
    num_list = [int(i) for i in numbers.split(delim)]
    if len(num_list) == 1:
        return(num_list[0])
    elif len(num_list) > 1:
        return(round(sum(num_list) / len(num_list)))


def report_summary_string(inlist, delim=" "):
    """Report summary string from a list of values

    E.g. 5 5 5 2 3 2 would be summarized as 5*3 2*2 3*1

    Parameters
    ----------
    inlist : list
        List of strings or values that can be coerced as strings, to be 
        summarized
    delim : str
        Character to separate the individual values
    """
    inlist = [str(i) for i in inlist]  # Convert to str
    counts = {i: inlist.count(i) for i in set(inlist)}
    # sort items by descending order of counts
    outlist = [str(i[0])+"*"+str(i[1])
               for i in sorted(counts.items(), reverse=True, key=lambda item: item[1])]
    return(delim.join(outlist))


def report_list_modes(inlist):
    """Report the mode of an input list

    Parameters
    ----------
    counts : list
        List of values that can be counted

    Returns
    -------
    list
        List of mode values. If there is a tie then the list length > 1
    """
    counts = {i: inlist.count(i) for i in set(inlist)}
    # check for ties
    maxvals = [i for i in counts if counts[i] == max(counts.values())]
    return(maxvals)


def get_not_gaps(start, end, gaplist):
    """Get inverse coords of gaps

    Parameters
    ----------
    start : int
    end : int
        Coords for original sequence, 0-based pythonic
    gaplist : list
        list of (int, int) tuples specifying gaps. Must not be overlapping
        otherwise GIGO

    Returns
    -------
    list
        list of (int, int) tuples specifying not-gaps
    """
    coords = [start]
    sortedgaps = sorted(gaplist, key=lambda x: int(x[0]))
    for tup in sortedgaps:
        coords.extend([tup[0], tup[1]])
    coords.append(end)
    # Reslice
    out = []
    for i in range(int(len(coords)/2)):
        out.append((coords[2*i], coords[2*i+1]))
    return(out)


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

    def __iter__(self):
        """Return iterator object for GFF.
        Iterate over keys of the _gffDict object
        """
        return(iter(self._gffDict.keys()))

    def __len__(self):
        """Length of a Gff object.
        Same as length of internal _gffDict object
        """
        return(len(self._gffDict))

    def addEntry(self, linearr, gffid):
        """ Add single GFF entry to Gff object

        Parameters
        ----------
        linearr : list
            List of GFF3 fields, length 9
        gffid : str
            ID attribute of the entry
        """
        # Check that GFF line has only nine fields
        if len(linearr) != 9:
            raise Exception(
                f'GFF3 input encountered with incorrect number of fields: {linearr}')
        idval = ""
        idsearch = re.search(r"ID=([^;]+)", linearr[8])
        if idsearch:
            idval = idsearch.group(1)
        elif gffid:
            idval = gffid
        else:
            raise Exception('GFF3 input without ID field')
        # Make a dict of each GFF3 column name and value
        self._gffDict[idval] = dict(zip(SharedValues.GFF3COLUMNS, linearr))
        # Split attributes field into key-value pairs and make a dict of it
        attrdict = {attr.split("=")[0] : attr.split("=")[1] for attr in linearr[8].rstrip(";").split(";")}
        self._gffDict[idval]['attrdict'] = attrdict

    def getValue(self, gffid, column):
        """Get value of column for a given GFF entry.
        GFF entry is specified by the ID.
        Returns the value.

        Parameters
        ----------
        gffid : str
            ID of the GFF entry
        column : str
            Name of the column to retrieve
        """
        if column in SharedValues.GFF3COLUMNS:
            if gffid in self._gffDict:
                return(self._gffDict[gffid][column])
            else:
                raise Exception("Unknown GFF3 ID " + gffid)
        else:
            raise Exception("Unknown GFF3 column name " + column)

    def changeValue(self, gffid, column, newvalue):
        """Change value of column for a given GFF entry.
        GFF entry is specified by the ID.

        Parameters
        ----------
        gffid : str
            ID of the GFF entry
        column : str
            Name of the column to retrieve
        newvalue
            New value of column to put in entry
        """
        if column in SharedValues.GFF3COLUMNS:
            if gffid in self._gffDict:
                self._gffDict[gffid][column] = newvalue
            else:
                raise Exception("Unknown GFF3 ID " + gffid)
        else:
            raise Exception("Unknown GFF3 column name " + column)

    def getAttr(self, gffid, attribute):
        """Get value from attributes field of a GFF entry.

        Parameters
        ----------
        gffid : str
            ID of the GFF entry
        attribute : str
            Key of the attribute requested
        """
        if gffid in self._gffDict:
            if attribute in self._gffDict[gffid]['attrdict']:
                return(self._gffDict[gffid]['attrdict'][attribute])
            else:
                logger.debug(
                    f"Attribute {attribute} not found for GFF ID {gffid}")
                return(None)
        else:
            raise Exception("Unknown GFF3 ID " + gffid)

    def changeAttr(self, gffid, attribute, newvalue):
        """Change value from attributes field of a GFF entry.

        Parameters
        ----------
        gffid : str
            ID of the GFF entry
        attribute : str
            Key of the attribute requested
        newvalue
            New value to put in the field
        """
        if gffid in self._gffDict:
            if attribute in self._gffDict[gffid]['attrdict']:
                self._gffDict[gffid]['attrdict'][attribute] = newvalue
                # Change attributes field
                attrlist = [ attr + "=" + str(self._gffDict[gffid]['attrdict'][attr])
                                for attr in self._gffDict[gffid]['attrdict']]
                self._gffDict[gffid]['attributes'] = ";".join(attrlist)
            else:
                logger.debug(
                    f"Unknown attribute {attribute} for GFF ID {gffid}")
        else:
            raise Exception("Unknown GFF3 ID " + gffid)

    def getEntry(self, gffid):
        """Get entry from Gff object with ID key, returns a list

        Parameters
        ----------
        gffid : str
            ID of the GFF entry
        """
        linearr = [self._gffDict[gffid][colname]
                   for colname in SharedValues.GFF3COLUMNS]
        return(linearr)

    def combineGff(self, gff):
        """Combine another Gff object into the current one by appending entries

        Checks for conflicting keys, raises exception if ID attribute in the
        second Gff object is already present in the current one.

        Parameters
        ----------
        gff : Gff
            Another GFF object to merge into the present one
        """
        for gffid in gff:
            if gffid not in self._gffDict:
                self.addEntry(gff.getEntry(gffid), gffid)
            else:
                raise Exception(
                    f"Attempting to merge two Gff objects with same ID {gffid} present in both")

    def list2gff(self, gfflist):
        """Add entries to Gff object from a list of GFF3 lines.

        Parameters
        ----------
        gfflist : list
            List of GFF3 records (list of str)
        """
        # Iterate through list, split each line into columns
        for line in gfflist:
            if not re.match(r"#", line):  # Skip header lines
                line = line.rstrip()  # Strip trailing whitespace
                linearr = line.split("\t")
                self.addEntry(linearr, None)

    def gff2list(self):
        """Write Gff3 object to list of strings, for printing.
        Returns list of str representing GFF3 entries, sorted.
        Effectively reconstitutes the input to list2gff()
        """

        outarr = []
        for gffid in self._gffDict:  # ID value
            linearr = self.getEntry(gffid)
            # Make numeric fields int - necessary for correct sorting later
            linearr[3] = int(linearr[3])
            linearr[4] = int(linearr[4])
            outarr.append(linearr)
        # Sort by seqid, start, and end fields in that order
        outarr = sorted(outarr, key=itemgetter(0, 3, 4))
        # Covnert fields back to str and join back into a tab-separated string
        outarr = [map(str, linearr) for linearr in outarr]
        outarr = ["\t".join(linearr) for linearr in outarr]
        return(outarr)

    def file2gff(self, filename):
        """Read in GFF3 file directly to Gff object

        Parameters
        ----------
        filename : str
            Path to file to open
        """
        with open(filename, "r") as fh:
            slurp = [line.rstrip() for line in fh]
            self.list2gff(slurp)

    def gff2file(self, filename, header=True):
        """Write Gff object directly to GFF3 file

        Parameters
        ----------
        filename : str
            Path to file to write
        header : bool
            Include ##gff-version header?
        """
        fh = open(filename, "w")
        outarr = self.gff2list()
        if header:
            fh.write("##gff-version 3\n")
        for line in outarr:
            fh.write(line+"\n")
        fh.close()

    def gff2fh(self, fh, header=True):
        """Write Gff object to open filehandle

        Parameters
        ----------
        fh : _io.TextIOWrapper
            File handle
        header : bool
            Include ##gff-version header?
        """
        outarr = self.gff2list()
        if header:
            fh.write("##gff-version 3\n")
        for line in outarr:
            fh.write(line+"\n")
