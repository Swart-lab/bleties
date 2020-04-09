#!/usr/bin/env python3

import re
import sys
import logging
import pysam
import json

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
    # Report putative IESs as list of GFF records and list of SeqRecord objects
    logging.info("Reporting putative IESs in GFF format")
    (iesgff, iesseq) = iesrecords.reportPutativeIes(args.min_break_coverage, args.min_del_coverage)
    # Write gff version header and command line as comment
    args.out.write("##gff-version 3\n")
    args.out.write("# " + " ".join(sys.argv) + "\n")
    # Write each GFF entry as a tab-separated line
    iesgff.gff2fh(args.out)
    # Write Fasta file of putative IES sequences
    if args.out_fasta:
        logging.info("Reporting consensus sequences of putative IESs to Fasta file "+args.out_fasta)
        SeqIO.write(iesseq, args.out_fasta, "fasta")
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
