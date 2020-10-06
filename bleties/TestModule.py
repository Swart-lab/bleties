#!/usr/bin/env python3

import unittest

from bleties import Milraa
from bleties import SharedFunctions
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
        cigartuples = [(4,2), (0,3), (1,3), (2,3), (4,3)]
        rstart = 10
        self.assertEqual(
                SharedFunctions.getCigarOpQuerySeqs(
                    qseq, cigartuples, rstart, target_op="S"),
                [("AA",0,2,10,10),("TTA",8,11,16,16)])

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
        seqlist = ['ATCG','ATCG','ATGG']
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


# class TestMilret(unittest.TestCase):


if __name__ == '__main__':
    unittest.main()
