#!/usr/bin/env python3

import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

from bleties.SharedFunctions import Gff


logger = logging.getLogger("Insert")


class Insert(object):
    """Insert IESs from GFF file and associated Fasta file to a reference
    assembly to produce MAC+IES genome sequence"""

    def __init__(self, refgenome, gff, iesfasta):
        """Initialize Insert object

        Parameters
        ----------
        refgenome : dict
            dict of SeqRecord objects, keyed by sequence ID, representing the
            MAC reference genome without IESs.
        gff : Gff
            GFF containing coordinates of IESs to be inserted into reference
        iesfasta : dict
            dict of SeqRecord objects, keyed by sequence ID, representing IES
            sequences to be inserted into the reference MAC genome.
        """
        self._refgenome = refgenome
        # Make copy of refgenome to modify, otherwise original will also be
        # modified
        self._newgenome = {ctg : SeqRecord(Seq(str(refgenome[ctg].seq)),
                                           id=ctg, name=refgenome[ctg].name,
                                           description=refgenome[ctg].description)
                            for ctg in refgenome}
        self._gff = gff
        self._iesfasta = iesfasta
        self._iesdict = defaultdict( # contig
                lambda: defaultdict( # pos
                    dict))           # dict of gffid, seqlen, seq, newstart, newend
        self._newgff = Gff()        # Gff entries for IESs with updated coords


    def _filterInserts(self):
        """Process GFF file and filter putative IESs to be inserted

        If multiple inserts at the same location, choose the longer one
        If sequence contains X or - characters, also reject
        """
        for gffid in self._gff:
            # Check for junctions only, i.e. start == end
            start = self._gff.getValue(gffid, 'start')
            end = self._gff.getValue(gffid, 'end')
            if start == end:
                ctg = self._gff.getValue(gffid, 'seqid')
                if gffid in self._iesfasta:
                    seq = str(self._iesfasta[gffid].seq)
                    # reject if there are X or - characters
                    if 'X' in seq or '-' in seq:
                        logger.debug(f"Skipping sequence {gffid} because of X and - characters")
                    else:
                        # if an insert already recorded, take the longer one
                        if self._iesdict[ctg][start]:
                            if len(seq) > self._iesdict[ctg][start]['seqlen']:
                                logger.debug(f"More than one insert at location {ctg} {str(start)}, choosing the longer one {gffid}")
                                self._iesdict[ctg][start] = {
                                        'seq' : seq,
                                        'seqlen' : len (seq),
                                        'gffid' : gffid,
                                        'newstart' : int(start),
                                        'newend' : int(end)}
                        else:
                            self._iesdict[ctg][start] = {
                                    'seq' : seq,
                                    'seqlen' : len(seq),
                                    'gffid' : gffid,
                                    'newstart' : int(start),
                                    'newend' : int(start) + len(seq)}
                else:
                    logger.warn(f"Insert sequence ID {gffid} not in Fasta file!")


    def _updatePositions(self):
        """After filtering inserts and recording them in iesdict, update their
        coordinates after adding the sequences in
        
        Run this after filterInserts()
        """
        for ctg in self._iesdict:
            sp = sorted(self._iesdict[ctg], key=lambda i: int(i))
            for i in range(len(sp)):
                seqlen = int(self._iesdict[ctg][sp[i]]['seqlen'])
                for j in range(i, len(sp)):
                    oldstart = int(self._iesdict[ctg][sp[j]]['newstart'])
                    oldend = int(self._iesdict[ctg][sp[j]]['newend'])
                    gffid = self._iesdict[ctg][sp[j]]['gffid']
                    newstart = oldstart + seqlen
                    newend = oldend + seqlen
                    # New coordinates in dict are 0-based [), python convention
                    self._iesdict[ctg][sp[j]]['newstart'] = newstart
                    self._iesdict[ctg][sp[j]]['newend'] = newend
                    # Record new Gff entries with updated columns
                    # start coordinate needs to be changed to follow GFF
                    # convention
                    self._newgff.addEntry(self._gff.getEntry(gffid), gffid)
                    self._newgff.changeValue(gffid, 'start', newstart + 1)
                    self._newgff.changeValue(gffid, 'end', newend)


    def _addSequences(self):
        """Insert IES sequences to the reference assembly
        
        Run this after updatePositions()
        """
        for ctg in self._iesdict:
            sp = sorted(self._iesdict[ctg], key=lambda i: int(i))
            for i in sp:
                insseq = self._iesdict[ctg][i]['seq']
                # Use coordinates from iesdict because these follow python
                # convention
                inspos = int(self._iesdict[ctg][i]['newstart'])
                self._newgenome[ctg].seq = self._newgenome[ctg].seq[0:inspos] + insseq + self._newgenome[ctg].seq[inspos:]


    def reportModifiedReference(self):
        """Add IESs to reference genome, report modified MAC+IES assembly and
        GFF file containing updated positions
        
        Returns
        -------
        dict
            dict of SeqRecords representing modified reference genome
        Gff
            IES records with updated inserted positions
        """
        self._filterInserts()
        self._updatePositions()
        ctg_orig_lengths = defaultdict(int)
        for ctg in self._refgenome:
            ctg_orig_lengths[ctg] = len(self._refgenome[ctg].seq)
        self._addSequences()
        ctg_new_lengths = defaultdict(int)
        for ctg in self._newgenome:
            ctg_new_lengths[ctg] = len(self._newgenome[ctg].seq)
        oldtotal = sum(ctg_new_lengths.values())
        newtotal = sum(ctg_orig_lengths.values())
        addedlen = newtotal - oldtotal
        logging.info(f"Original contigs total length: {str(oldtotal)}")
        logging.info(f"Modified contigs total length: {str(newtotal)}")
        logging.info(f"Total sequence length added: {str(addedlen)}")
        return(self._newgenome, self._newgff)
