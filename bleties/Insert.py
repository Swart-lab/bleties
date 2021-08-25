#!/usr/bin/env python3

import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

from bleties.SharedFunctions import Gff, get_not_gaps


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
        self._newgenome = {ctg: SeqRecord(Seq(str(refgenome[ctg].seq)),
                                          id=ctg, name=refgenome[ctg].name,
                                          description=refgenome[ctg].description)
                           for ctg in refgenome}
        self._gff = gff
        self._iesfasta = iesfasta
        self._iesdict = defaultdict(  # contig
            lambda: defaultdict(  # pos
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
                        logger.debug(
                            f"Skipping sequence {gffid} because of X and - characters")
                    else:
                        # if an insert already recorded, take the longer one
                        if self._iesdict[ctg][start]:
                            if len(seq) > self._iesdict[ctg][start]['seqlen']:
                                logger.debug(
                                    f"More than one insert at location {ctg} {str(start)}, choosing the longer one {gffid}")
                                self._iesdict[ctg][start] = {
                                    'seq': seq,
                                    'seqlen': len(seq),
                                    'gffid': gffid,
                                    'oldstart': int(start),
                                    'oldend': int(end) + len(seq),
                                    'newstart': int(start),
                                    'newend': int(end) + len(seq)
                                    }
                        else:
                            self._iesdict[ctg][start] = {
                                'seq': seq,
                                'seqlen': len(seq),
                                'gffid': gffid,
                                'oldstart': int(start),
                                'oldend': int(start) + len(seq),
                                'newstart': int(start),
                                'newend': int(start) + len(seq)
                                }
                else:
                    logger.warn(
                        f"Insert sequence ID {gffid} not in Fasta file!")

    def _filterDeletions(self):
        """Process GFF file and filter putative IESs to be deleted

        Overlapping regions will be skipped, only non-overlapping features
        retained

        Returns
        -------
        dict
            Dict keyed by contig ID (seqid), values are lists of dicts, each 
            containing details on a specific deletion, with keys seqid, start
            end, gffid. start and end are 0-based pythonic coordinates

        """
        coords = defaultdict(list)
        filtcoords = defaultdict(list)
        # Get region features only
        for gffid in self._gff:
            # Check for junctions only, i.e. start == end
            start = int(self._gff.getValue(gffid, 'start'))
            end = int(self._gff.getValue(gffid, 'end'))
            if start < end:
                seqid = self._gff.getValue(gffid, 'seqid')
                # Convert coords to 0-based python convention
                coords[seqid].append(
                    {'seqid': seqid, 'start': start-1, 'end': end, 'gffid': gffid})
        # Check for overlaps
        for seqid in coords:
            sortrecs = sorted(coords[seqid], key=lambda x: int(x['start']))
            currend = None
            for i in range(len(sortrecs)):
                gffid = sortrecs[i]['gffid']
                if currend:
                    if currend > int(sortrecs[i]['start']):
                        currend = max([currend, int(sortrecs[i]['end'])])
                        logger.debug(f"Overlapping record {gffid} skipped")
                    else:
                        currend = int(sortrecs[i]['end'])
                        if i < len(sortrecs) - 1:
                            if currend < int(sortrecs[i+1]['start']):
                                filtcoords[seqid].append(sortrecs[i])
                            else:
                                logger.debug(
                                    f"Overlapping record {gffid} skipped")
                        else:
                            # Sweep up last entry
                            filtcoords[seqid].append(sortrecs[i])
                else:
                    currend = int(sortrecs[i]['end'])
                    if i < len(sortrecs)-1:
                        if currend < int(sortrecs[i+1]['start']):
                            filtcoords[seqid].append(sortrecs[i])
                        else:
                            logger.debug(f"Overlapping record {gffid} skipped")
                    else:
                        # Sweet up last first entry
                        filtcoords[seqid].append(sortrecs[i])
        return(filtcoords)

    def _updatePositionsInserts(self):
        """After filtering inserts and recording them in iesdict, update their
        coordinates after adding the sequences in

        Run this after filterInserts()
        """
        for ctg in self._iesdict:
            sp = sorted(self._iesdict[ctg], key=lambda i: int(i))
            for i in range(len(sp)):
                seqlen = int(self._iesdict[ctg][sp[i]]['seqlen'])
                for j in range(i, len(sp)):
                    gffid = self._iesdict[ctg][sp[j]]['gffid']
                    if j > i:  # don't add if this is the same entry
                        # New coordinates in dict are 0-based [), python convention
                        self._iesdict[ctg][sp[j]]['newstart'] += seqlen
                        self._iesdict[ctg][sp[j]]['newend'] += seqlen
                    # Record new Gff entries with updated columns
                    # start coordinate needs to be changed to follow GFF
                    # convention
                    self._newgff.addEntry(self._gff.getEntry(gffid), gffid)
                    self._newgff.changeValue(gffid, 'start', self._iesdict[ctg][sp[j]]['newstart'] + 1)
                    self._newgff.changeValue(gffid, 'end', self._iesdict[ctg][sp[j]]['newend'])


    def _updatePointerPositionsInserts(self):
        """After updating coordinates of inserts, check for IES features that
        contain TA or pointer coordinates in attributes field, which may or may
        not differ from the start/end position, and update these too.
        
        Run this after _updatePositionsInserts()
        """
        for ctg in self._iesdict:
            for start in self._iesdict[ctg]:
                gffid = self._iesdict[ctg][start]['gffid']
                offset = int(self._iesdict[ctg][start]['newstart']) - int(self._iesdict[ctg][start]['oldstart'])
                # print(gffid) # testing
                # print(f"Offset {str(offset)}") # testing
                #
                # Check for TA or pointer coordinates in attributes of GFF and 
                # update with offset
                tps = self._newgff.getAttr(gffid, 'ta_pointer_start')
                pps = self._newgff.getAttr(gffid, 'pp_pointer_start')
                newend = self._iesdict[ctg][start]['newend']
                if tps:
                    # plus 1 to starts because 'newstart' and 'newend' are
                    # 0-based end exclusive coords, only converted to GFF
                    # convention in self._newgff
                    tps_offset = int(tps) - int(self._iesdict[ctg][start]['oldstart'])
                    # print(f"tps offset {str(tps_offset)}") # testing
                    # print(f"newstart {self._iesdict[ctg][start]['newstart']}")
                    # print(f"newend {self._iesdict[ctg][start]['newend']}")
                    self._newgff.changeAttr(gffid, 'ta_pointer_start', int(self._iesdict[ctg][start]['newstart']) + 1 + tps_offset)
                    self._newgff.changeAttr(gffid, 'ta_pointer_end', int(newend) + tps_offset)
                if pps:
                    pps_offset = int(pps) - int(self._iesdict[ctg][start]['oldstart'])
                    self._newgff.changeAttr(gffid, 'pp_pointer_start', int(self._iesdict[ctg][start]['newstart']) + 1 + pps_offset)
                    self._newgff.changeAttr(gffid, 'pp_pointer_end', int(newend) + pps_offset)


    def _updatePointerPositionsDeletions(self, dels):
        """After updating coordinates of deletions, check for IES features that
        contain TA or pointer coordinates in attributes field, which may or may
        not differ from the start/end position, and update these too.
        
        Run this after _updatePositionsDeletions()
        """
        for ctg in dels:
            for i in dels[ctg]:
                gffid = i['gffid']
                offset = int(i['newpos']) - int(i['start'])
                # Check for TA or pointer coordinates in attributes of GFF and 
                # update with offset
                tps = self._newgff.getAttr(gffid, 'ta_pointer_start')
                pps = self._newgff.getAttr(gffid, 'pp_pointer_start')
                if tps:
                    # minus 1 because tps is with GFF coords and oldstart on python coords
                    tps_offset = int(tps) - int(i['start']) - 1
                    # print(f"tps_offset {str(tps_offset)}") # testing
                    self._newgff.changeAttr(gffid, 'ta_pointer_start', int(i['newpos']) + tps_offset)
                    # when IES is deleted, TA start == TA end because this is
                    # now a zero-length junction-type feature
                    self._newgff.changeAttr(gffid, 'ta_pointer_end', int(i['newpos']) + tps_offset)
                if pps:
                    pps_offset = int(pps) - int(i['start']) - 1
                    self._newgff.changeAttr(gffid, 'pp_pointer_start', int(i['newpos']) + pps_offset)
                    self._newgff.changeAttr(gffid, 'pp_pointer_end', int(i['newpos']) + pps_offset)


    def _updatePositionsDeletions(self, dels):
        """After filtering deletions and recording them, update their
        coordinates in Gff object after removing the sequences from contigs

        Run this after filterDeletions()
        """
        for seqid in dels:
            sortedrecs = sorted(dels[seqid], key=lambda x: int(x['start']))
            for i in sortedrecs:
                i['newpos'] = i['start']
            for i in range(len(sortedrecs)):
                inslen = sortedrecs[i]['end'] - sortedrecs[i]['start']
                for j in range(i+1, len(sortedrecs)):
                    sortedrecs[j]['newpos'] = sortedrecs[j]['newpos'] - inslen
            # Record new Gff entry and update coordinates
            for i in sortedrecs:
                self._newgff.addEntry(
                    self._gff.getEntry(i['gffid']), i['gffid'])
                self._newgff.changeValue(i['gffid'], 'start', i['newpos'])
                self._newgff.changeValue(i['gffid'], 'end', i['newpos'])


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
                self._newgenome[ctg].seq = self._newgenome[ctg].seq[0:inspos] + \
                    insseq + self._newgenome[ctg].seq[inspos:]


    def _deleteSequences(self, dels):
        """
        Run this after filterDeletions()

        Parameters
        ----------
        dels : dict
            Output from _filterDeletions()
        """
        for seqid in dels:
            gaps = [(i['start'], i['end']) for i in dels[seqid]]
            end = len(self._refgenome[seqid].seq)
            notgaps = get_not_gaps(0, end, gaps)
            newseq = ""
            for i in notgaps:
                newseq += str(self._refgenome[seqid].seq[i[0]:i[1]])
            self._newgenome[seqid].seq = Seq(newseq)


    def reportInsertedReference(self):
        """Add IESs to reference genome, report modified MAC+IES assembly and
        GFF file containing updated positions

        Only IESs defined in the GFF as inserts (i.e. start == end) will be 
        added to the reference sequence. IES sequences containing gaps or
        mismatches (- or X characters) will be skipped.

        Returns
        -------
        dict
            dict of SeqRecords representing modified reference genome
        Gff
            IES records with updated region coords. Only the start and end
            columns are changed, and the attributes for adjusted coordinates
            relative to putative pointers or TA junctions: ta_pointer_start,
            ta_pointer_end, pp_pointer_start, pp_pointer_end.

            All other columns are inherited from the original GFF records
        """
        self._filterInserts()
        self._updatePositionsInserts()
        self._updatePointerPositionsInserts() # update ta_pointer... pp_pointer...
        self._addSequences()
        # Count difference in size
        oldtotal = sum([len(self._refgenome[ctg].seq)
                        for ctg in self._refgenome])
        newtotal = sum([len(self._newgenome[ctg].seq)
                        for ctg in self._newgenome])
        addedlen = newtotal - oldtotal
        logging.info(f"Original contigs total length: {str(oldtotal)}")
        logging.info(f"Modified contigs total length: {str(newtotal)}")
        logging.info(f"Total sequence length added: {str(addedlen)}")
        return(self._newgenome, self._newgff)


    def updateFeatureGff(self, feats):
        """Update coordinates of other annotations after filtering IES inserts

        Run this after _updatePositionsInserts().This function only works with
        insert mode, because in delete mode, the truncation or excision of
        features with potential annotated functions must be curated manually by
        the user.

        Parameters
        ----------
        feats : list
            list of lists; each element represents single GFF3 line split on
            tabs into columns

        Returns
        -------
        list
            list of lists; each element represents single GFF3 line. Do not use
            the SharedFunctions.Gff class because that implementation does not
            account for multi-segment features with the same ID but split
            across several lines.
        """
        newgff = []
        for ctg in self._iesdict:
            for feat in [f for f in feats if f[0] == ctg]:
                # Get ID from attributes
                attrdict = {attr.strip().split('=')[0] : attr.strip().split('=')[1]
                            for attr in feat[8].split(';') if attr}
                featid = attrdict['ID']
                # add up all offsets before this feature
                offset = sum([int(self._iesdict[ctg][i]['seqlen']) 
                              for i in self._iesdict[ctg] 
                              if int(i) < int(feat[3])])
                # check if any IES inserts fall inside the feature
                istart = int(feat[3])
                iend = int(feat[4])
                # get start positions of each IES feature that is inside
                # this feature
                iesinsides = [int(i) for i in self._iesdict[ctg]
                              if int(i) >= istart and int(i) < iend]
                iesinsides = sorted(iesinsides)
                insidelen = sum([int(self._iesdict[ctg][str(i)]['seqlen'])
                                 for i in iesinsides])
                # new start position
                newints = [istart + offset]
                # update coordinates with complement of new IES feature
                # positions
                for i in iesinsides:  # already sorted
                    # the following are already in GFF-coordinates
                    newints.append(self._iesdict[ctg][str(i)]['newstart'])
                    newints.append(self._iesdict[ctg][str(i)]['newend'] + 1)
                # new end position
                newints.append(iend + insidelen + offset)
                # split to tuples of start-end coords for non-IES complement
                # ranges
                newints_tuples = []
                for i in range(int(len(newints)/2)):
                    newints_tuples.append((newints[i*2], newints[i*2 + 1]))
                if len(newints_tuples) > 1:
                    # feature is interrupted by IESs, split into segments
                    for j in range(len(newints_tuples)):
                        newfeatid = featid + '.seg_' + str(j)
                        # Updates ID field if present, adds it if not
                        attrdict['ID'] = newfeatid
                        newattr = ';'.join([str(key) + '=' + str(attrdict[key])
                                            for key in attrdict])
                        # new GFF entry as list
                        newfeat = feat[0:3] + \
                            [newints_tuples[j][0], newints_tuples[j][1]] + \
                            feat[5:8] + [newattr]
                        # add entry to new Gff object
                        newgff.append(newfeat)
                else:
                    # new interval comprises only the new start and end, 
                    # the feature is not broken up into segments
                    # feature ID remains the same
                    newfeat = feat[0:3] + \
                        [newints_tuples[0][0], newints_tuples[0][1]] + \
                        feat[5:9]
                    newgff.append(newfeat)
        return(newgff)

    def reportDeletedReference(self):
        """Remove IESs from reference MAC+IES genome, report modified MAC-IES
        assembly and GFF file containing updated positions

        Only IESs defined in the GFF as retentions (i.e. end > start) will be
        removed from the reference sequence. IESs that overlap with each other
        will be skipped: the input GFF should be manually curated to select
        which region should actually be excised.

        Returns
        -------
        dict
            dict of SeqRecords representing modified reference genome
        Gff
            IES records with updated insert positions. Only the start and end
            columns are changed, all other columns are inherited from the 
            original GFF records
        """
        dels = self._filterDeletions()
        self._updatePositionsDeletions(dels)
        self._updatePointerPositionsDeletions(dels)
        self._deleteSequences(dels)
        # Count difference in size
        oldtotal = sum([len(self._refgenome[ctg].seq)
                        for ctg in self._refgenome])
        newtotal = sum([len(self._newgenome[ctg].seq)
                        for ctg in self._newgenome])
        deletedlen = oldtotal - newtotal
        logging.info(f"Original contigs total length: {str(oldtotal)}")
        logging.info(f"Modified contigs total length: {str(newtotal)}")
        logging.info(f"Total sequence length removed: {str(deletedlen)}")
        return(self._newgenome, self._newgff)
