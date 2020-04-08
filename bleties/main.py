#!/usr/bin/env python3

import re
import sys
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bleties import *

def milraa(args):
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
    alnfile = pysam.AlignmentFile(aln_filename, aln_mode)
    # Read reference Fasta file into memory
    refgenome = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
    # Initialize new IesRecords object to store putative IESs
    iesrecords = IesRecords.IesRecords(alnfile, aln_format, refgenome)
    # Process alignment to find putative IESs 
    iesrecords.findPutativeIes(args.min_ies_length)

    if args.dump:
        sys.stderr.write(str(iesrecords) + "\n") # Print summary of IesRecords object
        print(iesrecords.dump()) # Dump data to check

    # Report putative IESs as list of GFF records and list of SeqRecord objects
    (iesgff, iesseq) = iesrecords.reportPutativeIes(args.min_break_coverage, args.min_del_coverage)

    # Write gff version header and command line as comment
    args.out.write("##gff-version 3\n")
    args.out.write("# " + " ".join(sys.argv) + "\n")
    # Write each GFF entry as a tab-separated line
    for rec in iesgff:
        args.out.write("\t".join(rec)+"\n")

    # Write Fasta file of putative IES sequences
    if args.out_fasta:
        SeqIO.write(iesseq, args.out_fasta, "fasta")

    # Close AlignmentFile
    alnfile.close()


