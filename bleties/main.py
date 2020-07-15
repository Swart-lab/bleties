#!/usr/bin/env python3

import re
import sys
import logging
import pysam
import json
import statistics as stats
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bleties import *

def milraa(args):
    logging.info("Started MILRAA")
    # Check that only either SAM or BAM specified
    aln_filename = "-"
    aln_format = "sam"
    aln_mode = "r"
    if args.sam:
        if args.bam:
            sys.exit("Error: Specify either SAM or BAM input, not both")
        else:
            aln_filename = args.sam
            aln_format = "sam"
            aln_mode = "r"
    if args.bam:
        if args.sam:
            sys.exit ("Error: Specify either SAM or BAM input, not both")
        else:
            aln_filename=args.bam
            aln_format = "bam"
            aln_mode = "rb"
    # Open SAM or BAM file 
    logging.info("Opening alignment file "+aln_filename)
    alnfile = pysam.AlignmentFile(aln_filename, aln_mode)
    logging.info("Alignment file contains " 
                 + str(alnfile.mapped) 
                 + " reads mapped to " 
                 + str(alnfile.nreferences)
                 + " reference sequences")
    # Read reference Fasta file into memory
    logging.info("Reading mapping sequence file to memory "+args.ref)
    refgenome = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
    # Initialize new IesRecords object to store putative IESs
    iesrecords = Milraa.IesRecords(alnfile, aln_format, refgenome)
    # Process alignment to find putative IESs 
    logging.info("Processing alignment to find putative IESs")
    iesrecords.findPutativeIes(args.min_ies_length)
    if args.dump:
        logging.info("Dumping data in JSON format to STDOUT")
        # sys.stderr.write(str(iesrecords) + "\n") # Print summary of IesRecords object
        print(iesrecords.dump()) # Dump data to check
    # Report putative IESs as list of GFF records and dict of SeqRecord objects
    logging.info("Reporting putative IESs in GFF format")
    (iesgff, iesseq) = iesrecords.reportPutativeIes(args.min_break_coverage, args.min_del_coverage)

    # TEST: Report mismatch percentages of reads with and without each putative IES
    logging.info("Reporting possibly spurious IESs due to misassembly or mapped paralogs")
    test_type = 'mann-whitney' # TODO: user option
    # test_type = 't'
    with open("mismatch_pc.tsv","w") as fh_mm:
        fh_mm.write("\t".join(['ID',
            'mean_mismatch_reads_with_indel',
            'mean_mismatch_without_indel',
            'stdev_with_indel',
            'stdev_without_indel',
            'statistic',
            'p-value',
            'no_reads_with_indel',
            'no_reads_without_indel',
            'diagnosis']))
        fh_mm.write("\n")
        for bpid in iesgff:
            ins_mm, non_mm = iesrecords.reportIndelReadMismatchPc(
                iesgff.getValue(bpid,'seqid'),
                int(iesgff.getValue(bpid,'start')),
                int(iesgff.getValue(bpid,'end')),
                int(iesgff.getAttr(bpid,'IES_length'))
            )
            if test_type == 'mann-whitney':
                # Mann-Whitney U test for whether mismatch % with indel of interest
                # is greater than without
                mwstat, mwpval = mannwhitneyu(ins_mm, non_mm, alternative='greater')
            else:
                # Ward's t-test (non-equal population variances)
                mwstat, mwpval = ttest_ind(ins_mm, non_mm, equal_var=False)
            # Report
            outarr = [bpid, 
                round(stats.mean(ins_mm),2),
                round(stats.mean(non_mm),2),
                round(stats.stdev(ins_mm),2),
                round(stats.stdev(non_mm),2),
                round(mwstat,2),
                '%.2E' % mwpval, # scientific notation
                len(ins_mm),
                len(non_mm)]
            diagnosis = "ok"
            if len(ins_mm) > len(non_mm):
                diagnosis = 'misassembly?'
            PVAL_UNCORR = 0.05 # TODO magic number
            pval_corr = PVAL_UNCORR / len(iesgff) # Bonferroni correction
            if mwpval < pval_corr and stats.mean(ins_mm) > stats.mean(non_mm):
                diagnosis = "paralog?"
            outarr.append(diagnosis)
            fh_mm.write("\t".join([str(i) for i in outarr]))
            fh_mm.write("\n")
            # fh_mm.write(" ".join([str(i) for i in ins_mm]) + "\n")
            # fh_mm.write(" ".join([str(i) for i in non_mm]) + "\n")

    # Write gff version header and command line as comment
    args.out.write("##gff-version 3\n")
    args.out.write("# " + " ".join(sys.argv) + "\n")
    # Write each GFF entry as a tab-separated line
    iesgff.gff2fh(args.out)
    # Write Fasta file of putative IES sequences
    if args.out_fasta:
        logging.info("Reporting consensus sequences of putative IESs to Fasta file "+args.out_fasta)
        SeqIO.write(iesseq.values(), args.out_fasta, "fasta")
    # Report junction sequences
    if args.out_junction and args.junction_flank:
        junctionseqs = Milraa.getIndelJunctionSeqs(iesgff, iesseq, refgenome, args.junction_flank)
        logging.info("Reporting flanking sequences of putative IESs to file "+args.out_junction)
        with open(args.out_junction, "w") as fhjunc:
            fhjunc.write("\t".join(['id','leftflank','rightflank','indel','ref'])+"\n") # header
            for junc in junctionseqs:
                fhjunc.write("\t".join(junc) + "\n")
    # Close AlignmentFile
    alnfile.close()
    logging.info("Finished MILRAA")

def milret(args):
    logging.info("Started MILRET")
    # Read BAM file - SAM not supported because we need random access
    logging.info("Opening alignment file +"+args.bam)
    alnfile = pysam.AlignmentFile(args.bam, "rb") 
    # Initialize IesRetentionsMacOnly object
    iesretentions = Milret.IesRetentionsMacOnly(args.ies, alnfile)
    # Count mapping operations per site
    logging.info("Counting IES+ and IES- forms at each junction in file "+args.ies)
    iesretentions.findMappingOps()
    # Report retention scores to file
    logging.info("Calculating retention scores per junction")
    iesretentions.calculateRetentionScores()
    iesretentions.reportRetentionScores(args.out)
    logging.info("Finished MILRET")
