

import unittest
from bleties import Milret
from bleties import Milraa


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
        # Test region IESs
        self.assertEqual(
                Milraa.getPointers(seq, 2, 13, iesseq_TG, "test_del_left"),
                ("TAG", 2, 13))
        self.assertEqual(
                Milraa.getPointers(seq, 5, 16, iesseq_TG, "test_del_right"),
                ("TAG", 2, 13))


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


class TestMilret(unittest.TestCase):

    def test_getOperationAtRefPos(self):
        self.assertEqual(
                Milret.getOperationAtRefPos(17, 15, "20M77D12I77M", 1, 1),"M")
        self.assertEqual(
                Milret.getOperationAtRefPos(37, 15, "20M77D12I77M", 1, 1),"D")
        self.assertEqual(
                Milret.getOperationAtRefPos(97, 1, "20M77D12I77M", 1, 1),"I")
        self.assertEqual(
                Milret.getOperationAtRefPos(98, 1, "20M77D12I77M", 1, 1),"M")


if __name__ == '__main__':
    unittest.main()
