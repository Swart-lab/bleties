#!/usr/bin/env python3

import logging

from bleties.SharedFunctions import *
from bleties.SharedValues import SharedValues


logger = logging.getLogger("Milcor")


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
        self._perRead = defaultdict(  # read ID
            lambda: defaultdict(
                list))  # list of IES IDs
        # Initialize dict to record IES co-occurrences
        self._perIes = defaultdict(  # IES ID
            lambda: defaultdict(  # ID of other IESs that co-occur with first
                int))  # Number of co-occurrences/co-observations

    def countIesCooccurrences(self, match_lengths=False, threshold=0.05,
                              fetch_ctg=None, fetch_start=None, fetch_stop=None):
        """Iterate through alignment and count IES co-occurrences per read

        Parameters
        ----------
        match_lengths : bool
            Match lengths of detected inserts in reads with the reported IES
            length in GFF file?
        threshold : float
            Relative length difference to accept when matching IES lengths,
            inclusive. For example, if reported IES length is 100 bp and
            threshold is 0.05, accept insert lengths between 95 and 105
            inclusive as IES present. Only used if match_lengths is True
        fetch_ctg : str
            Name of contig to fetch alignments for. If None, all alignments used
        fetch_start : int
        fetch_stop : int
            Coordinates of contigs to fetch alignments for. If None, all
            alignments used. 1-based GFF convention
        """
        # Iterate across all reads
        if fetch_start:
            fetch_start = fetch_start - 1
        for alnrec in self._alnfile.fetch(contig=fetch_ctg, start=fetch_start, stop=fetch_stop):
            qname = alnrec.query_name
            if (not alnrec.is_secondary) and (not alnrec.is_supplementary):
                ctg = alnrec.reference_name
                start = int(alnrec.reference_start) + \
                    1  # convert to 1-based inclusive
                end = int(alnrec.reference_end)  # 1-based inclusive
                if not qname:
                    logger.warn(
                        "Query in alignment without a query_name field")
                if self._perRead[qname]['ref'] and self._perRead[qname]['ref'] != ctg:
                    logger.warn(f"Read {qname} aligns to more than one contig")
                # Record info on each aligned read
                self._perRead[qname]['ref'] = ctg
                self._perRead[qname]['start'] = start
                self._perRead[qname]['end'] = end
                # Check if this contig has IESs defined
                if ctg in self._gffByContig:
                    for rec in nested_dict_to_list_fixed_depth(self._gffByContig[ctg], 2):
                        iesstart, iesend, iesid_list = rec
                        # only count junctions, i.e. insertions relative to reference
                        if iesstart >= start and iesstart <= end and iesstart == iesend:
                            res = getOperationAtRefPos(iesstart,
                                                       start,
                                                       alnrec.cigarstring,
                                                       1, 1)
                            if res:
                                op, op_len = res
                                if op == "I":
                                    # Read has insert corresponding to IES position
                                    for iesid in iesid_list:
                                        if match_lengths:
                                            # Get IES length from GFF input
                                            # take average if more than one length present
                                            ieslength = self._gff.getAttr(
                                                iesid, "IES_length")
                                            if ieslength:  # None returned if attribute absent
                                                ieslength = mean_of_number_list(
                                                    ieslength)
                                                # Only add if length of "I" operation matches reported IES length
                                                if op_len >= ieslength*(1-threshold) and op_len <= ieslength*(1+threshold):
                                                    # print(f"Length of IES {str(op_len)} matches insert length {ieslength} on read {qname}")
                                                    self._perRead[qname]['present'].append(
                                                        iesid)
                                                else:
                                                    self._perRead[qname]['absent'].append(
                                                        iesid)
                                            else:
                                                logger.warn(
                                                    f"Input GFF file does not report IES length reported for IES {iesid}")
                                        else:
                                            # Ignore IES lengths, count all inserts at the position
                                            self._perRead[qname]['present'].append(
                                                iesid)
                                    # TODO also work for deletions
                                elif op == "M":
                                    # No IES in the defined position
                                    for iesid in iesid_list:
                                        self._perRead[qname]['absent'].append(
                                            iesid)

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
            fh.write(json.dumps({"_perRead": self._perRead,
                                 "_perIes": self._perIes}, indent=4))

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

    def binReads(self, macreads, micreads, otherreads, noiesreads, threshold=0.9,
                 fetch_ctg=None, fetch_start=None, fetch_stop=None):
        """Bin reads into MAC or MIC origin

        Parameters
        ----------
        macreads : str
            Path to write MAC reads in Fasta format
        micreads : str
            Path to write MIC reads in Fasta format
        otherreads : str
            Path to write ambiguous reads which do not cross threshold
        noiesreads : str
            Path to write reads that do not span IES sites
        threshold : float
            Minimum proportion of reads with IES excision/retention to classify
            as MAC or MIC respectively.
        """

        if threshold > 1 or threshold < 0:
            raise Exception("Binning threshold must be between 0 and 1")

        fh_mac = open(macreads, "w")
        fh_mic = open(micreads, "w")
        fh_oth = open(otherreads, "w")
        fh_non = open(noiesreads, "w")

        counter = 0
        if fetch_start:
            fetch_start = fetch_start - 1
        for rec in self._alnfile.fetch(contig=fetch_ctg, start=fetch_start, stop=fetch_stop):
            if (not rec.is_unmapped) and (not rec.is_secondary) and (not rec.is_supplementary):
                qname = rec.query_name
                qseq = rec.query_sequence
                if qname and qname in self._perRead:
                    counter += 1
                    if counter % 1000 == 0:
                        logger.debug(f"Processed {counter} reads")

                    iesplus = len(self._perRead[qname]['present'])
                    iesminus = len(self._perRead[qname]['absent'])
                    iestotal = iesplus + iesminus

                    if iestotal == 0:
                        # Non IES read
                        fh_non.write(f">{qname}\n")
                        fh_non.write(f"{qseq}\n")
                    else:
                        iesret = iesplus / iestotal
                        if iesret >= threshold:  # MIC read
                            # MIC read
                            fh_mic.write(f">{qname}\n")
                            fh_mic.write(f"{qseq}\n")
                        elif iesret <= 1-threshold:
                            # MAC read
                            fh_mac.write(f">{qname}\n")
                            fh_mac.write(f"{qseq}\n")
                        else:
                            # other read
                            fh_oth.write(f">{qname}\n")
                            fh_oth.write(f"{qseq}\n")
                else:
                    logger.debug(f"Read {qname} not in _perRead dict")

        fh_mac.close()
        fh_mic.close()
        fh_oth.close()
        fh_non.close()

    # TODO: ? Graphical representation of IES co-occurrence
    # TODO: Flag IESs that have co-occurrences more or less than expected
