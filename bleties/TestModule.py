#!/usr/bin/env python3

import unittest
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo

from bleties import Milraa
from bleties import SharedFunctions
from bleties import Insert
# from bleties import Milret

# To run tests from parent directory:
# python -m unittest -v bleties.TestModule


class TestSharedFunctions(unittest.TestCase):

    def test_getOperationAtRefPos(self):
        self.assertEqual(
            SharedFunctions.getOperationAtRefPos(17, 15, "20M77D12I77M", 1, 1),
            ("M", 20))
        self.assertEqual(
            SharedFunctions.getOperationAtRefPos(37, 15, "20M77D12I77M", 1, 1),
            ("D", 77))
        self.assertEqual(
            SharedFunctions.getOperationAtRefPos(97, 1, "20M77D12I77M", 1, 1),
            ("I", 12))
        self.assertEqual(
            SharedFunctions.getOperationAtRefPos(98, 1, "20M77D12I77M", 1, 1),
            ("M", 77))

    def test_mean_of_number_list(self):
        self.assertEqual(
            SharedFunctions.mean_of_number_list("12_13_14", "_"),
            13)

    def test_getCigarOpQuerySeqs(self):
        qseq = "AATACCCATTA"
        cigartuples = [(4, 2), (0, 3), (1, 3), (2, 3), (4, 3)]
        rstart = 10
        self.assertEqual(
            SharedFunctions.getCigarOpQuerySeqs(
                qseq, cigartuples, rstart, target_op="S"),
            [("AA", 0, 2, 10, 10), ("TTA", 8, 11, 16, 16)])

    def test_report_summary_string(self):
        inlist = [3, 5, 5, 5, 2, 2]
        self.assertEqual(
            SharedFunctions.report_summary_string(inlist, " "),
            "5*3 2*2 3*1")

    def test_report_list_modes(self):
        inlist = [1, 2, 2, 2, 2, 3, 4]
        self.assertEqual(
            SharedFunctions.report_list_modes(inlist),
            [2])
        inlist2 = ['+', '+', '+', '-']
        self.assertEqual(
            SharedFunctions.report_list_modes(inlist2),
            ['+'])
        inlist_tie = [1, 1, 2, 2]
        self.assertEqual(
            SharedFunctions.report_list_modes(inlist_tie),
            [1, 2])

    def test_get_not_gaps(self):
        start = 0
        end = 1000
        gaps = [(47, 51), (400, 640)]
        self.assertEqual(
            SharedFunctions.get_not_gaps(start, end, gaps),
            [(0, 47), (51, 400), (640, 1000)])


class TestMilraa(unittest.TestCase):

    def test_getPointers(self):
        seq = "ATAGCGCTGCGTTTAGTT"
        iesseq_tie = "TGAATGC"
        iesseq_TG = "TGAATAA"
        # Test insert IESs
        self.assertEqual(
            Milraa.getPointers(seq, 7, 7, iesseq_tie, "test_ins_tie"),
            ("tie", 7, 7))
        self.assertEqual(
            Milraa.getPointers(seq, 7, 7, iesseq_TG, "test_ins_TG"),
            ("TG", 7, 7))
        # Test deletion IESs
        self.assertEqual(
            Milraa.getPointers(seq, 2, 13, iesseq_TG, "test_del_left"),
            ("TAG", 2, 13))
        self.assertEqual(
            Milraa.getPointers(seq, 5, 16, iesseq_TG, "test_del_right"),
            ("TAG", 2, 13))

    def test_adjustPointerTA(self):
        # Test deletion IESs
        self.assertEqual(
            Milraa.adjustPointerTA(2, 13, "TTAG"),
            (3, 14, "TAG"))

    def test_adjustPointerMaxlength(self):
        seq = "ATAGCGCTGCGTTTAGTT"
        # Test insertion IESs
        ies = "AGCGT"
        self.assertEqual(
            Milraa.adjustPointerMaxlength(seq, 2, 2, "AG", ies),
            (1, 1, "TAG"))
        # Test deletion IESs
        self.assertEqual(
            Milraa.adjustPointerMaxlength(seq, 3, 14, "AG", None),
            (2, 13, "TAG"))

    def test_alnFromSeqs(self):
        seqlist = ['ATGCG',
                   'ATCG',
                   'ATCG']
        self.assertEqual(
            str(Milraa.alnFromSeqs(seqlist, 0.7).seq),
            "ATXCG")
        # Milraa.alnFromSeqs produces SeqRecord object, but direct comparison
        # of SeqRecord objects is deprecated, so must first convert to str

    def test_alnDumbFromSeqs(self):
        seqlist = ['ATCG', 'ATCG', 'ATGG']
        self.assertEqual(
            str(Milraa.alnDumbFromSeqs(seqlist, 0.7).seq),
            "ATXG")

    def test_getIndels(self):
        cigar = "5M2D3M2I5M"
        qseq = "ATATATTTCCATATA"
        self.assertEqual(
            Milraa.getIndels(cigar, 5, 2, qseq),
            [(10, 11, 0, "D", ""),
             (14, 14, 2, "I", "CC")])

    def test_alnRemoveGapOnlyCols(self):
        s1 = SeqRecord(Seq('A-TT---TTAA---'),id='s1',name='s1')
        s2 = SeqRecord(Seq('AATT---TTAA---'),id='s2',name='s2')
        aln = MultipleSeqAlignment([s1, s2])
        s1_nogap = SeqRecord(Seq('A-TTTTAA'),id='s1',name='s1')
        s2_nogap = SeqRecord(Seq('AATTTTAA'),id='s2',name='s2')
        alnnogap = MultipleSeqAlignment([s1_nogap, s2_nogap])
        aln = MultipleSeqAlignment([s1, s2])
        # Use format() to report, because the Align objects will be compared
        # by hash values which will not be equal
        self.assertEqual(
            Milraa.alnRemoveGapOnlyCols(aln).format('fasta'),
            alnnogap.format('fasta'))


# class TestMilret(unittest.TestCase):


class TestInsert(unittest.TestCase):
    # Class variables used by several tests
    ref = {'ctg1': SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAA'), id='ctg1'),
           'ctg2': SeqRecord(Seq('GGGGGGGGGGGGGGGGGGGG'), id='ctg2'),
           'ctg3': SeqRecord(Seq('CCCCCCCCCCTACCCCCCCC'), id='ctg3')
           }
    ies = {'ies1': SeqRecord(Seq('TTTT'), id='ies1'),
           'ies2': SeqRecord(Seq('GGGGG'), id='ies2'),
           'ies3': SeqRecord(Seq('CCCCC'), id='ies3'),
           'ies5': SeqRecord(Seq('TTTTT'), id='ies5'),
           'ies6': SeqRecord(Seq('CTAGG'), id='ies6')}
    gfflist = ["ctg1\t.\t.\t5\t5\t.\t.\t.\tID=ies1;",
               "ctg1\t.\t.\t9\t9\t.\t.\t.\tID=ies2;",
               "ctg2\t.\t.\t9\t9\t.\t.\t.\tID=ies3;",
               "ctg2\t.\t.\t15\t18\t.\t.\t.\tID=ies4;",
               "ctg3\t.\t.\t3\t3\t.\t.\t.\tID=ies5;",
               "ctg3\t.\t.\t9\t9\t.\t.\t.\tID=ies6;ta_pointer_start=11;ta_pointer_end=11;"
               ]
    oldfeatures = ["ctg1\t.\tgene\t3\t7\t.\t.\t.\tID=gene1;key1=attr1;key2=attr2",
                   "ctg1\t.\tgene\t12\t15\t.\t.\t.\tID=gene2"
                   ]


    def test_reportInsertedReference(self):
        gff = SharedFunctions.Gff()
        gff.list2gff(TestInsert.gfflist)
        ins = Insert.Insert(TestInsert.ref, gff, TestInsert.ies)
        newfasta, newgff = ins.reportInsertedReference()
        self.assertEqual(
            str(newfasta['ctg1'].seq),
            'AAAAATTTTAAAAGGGGGAAAAAAAAAAA')
        self.assertEqual(
            str(newfasta['ctg2'].seq),
            'GGGGGGGGGCCCCCGGGGGGGGGGG')


    def test_updateFeatureGff(self):
        iesgff = SharedFunctions.Gff()
        iesgff.list2gff(TestInsert.gfflist)
        ins = Insert.Insert(TestInsert.ref, iesgff, TestInsert.ies)
        ins._filterInserts()
        ins._updatePositionsInserts()
        annotgff = SharedFunctions.Gff()
        annotgff.list2gff(TestInsert.oldfeatures)
        newgff = ins.updateFeatureGff(annotgff)
        self.assertEqual(
            [str(i) for i in newgff.getEntry('gene1.seg_0')],
            ['ctg1','.','gene','3','5','.','.','.','ID=gene1.seg_0;key1=attr1;key2=attr2'])
        self.assertEqual(
            [str(i) for i in newgff.getEntry('gene1.seg_1')],
            ['ctg1','.','gene','10','11','.','.','.','ID=gene1.seg_1;key1=attr1;key2=attr2'])
        self.assertEqual(
            [str(i) for i in newgff.getEntry('gene2')],
            ['ctg1','.','gene','21','24','.','.','.','ID=gene2'])


    def test_updateFeatureGff_tapointer(self):
        iesgff = SharedFunctions.Gff()
        iesgff.list2gff(TestInsert.gfflist)
        ins = Insert.Insert(TestInsert.ref, iesgff, TestInsert.ies)
        ins._filterInserts()
        ins._updatePositionsInserts()
        ins._updatePointerPositionsInserts()
        ins._addSequences()
        self.assertEqual(str(ins._newgff.getValue('ies6', 'start')), '15')
        self.assertEqual(str(ins._newgff.getValue('ies6', 'end')), '19')
        self.assertEqual(str(ins._newgff.getAttr('ies6', 'ta_pointer_start')), '16')
        self.assertEqual(str(ins._newgff.getAttr('ies6', 'ta_pointer_end')), '20')


    def test_reportDeletedReference(self):
        gff = SharedFunctions.Gff()
        gff.list2gff(TestInsert.gfflist)
        # Insert sequences into reference
        ins = Insert.Insert(TestInsert.ref, gff, TestInsert.ies)
        newfasta, newgff = ins.reportInsertedReference()
        # print("\n".join(newgff.gff2list())) # testing
        # Take them out again
        dels = Insert.Insert(newfasta, newgff, None)
        delfasta, delgff = dels.reportDeletedReference()
        # print("\n".join(delgff.gff2list())) # testing
        # Check that they are the same sequence
        self.assertEqual(
            str(delfasta['ctg1'].seq),
            str(TestInsert.ref['ctg1'].seq))
        self.assertEqual(
            str(delfasta['ctg2'].seq),
            str(TestInsert.ref['ctg2'].seq))
        self.assertEqual(
            str(delfasta['ctg3'].seq),
            str(TestInsert.ref['ctg3'].seq))
        self.assertEqual(str(delgff.getAttr('ies6', 'ta_pointer_start')), '11')
        self.assertEqual(str(delgff.getAttr('ies6', 'ta_pointer_end')), '11')


if __name__ == '__main__':
    unittest.main()
