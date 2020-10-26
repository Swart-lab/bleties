#!/usr/bin/env python3

import re
import sys
import logging
import pysam
import json
import statistics as stats
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bleties import Milraa, Milret, Milcor, Miltel, Insert
from bleties import __version__
from bleties import SharedFunctions


def read_bam_ref(bam, fasta=None):
    """Read input SAM or BAM alignment and reference Fasta file.

    Checks for BAM index, and that headers in sequence files match those in BAM

    Parameters
    ----------
    bam : str
        Path to BAM file mapping
    fasta : str
        Path to Fasta file containing reference used for mapping, optional.

    Returns
    -------
    pysam.AlignmentFile
        Handle to BAM mapping file
    dict
        dict of Bio.SeqRecords objects representing reference genome
    """
    logger = logging.getLogger("main.read_bam_ref")

    # Open BAM file
    logger.info(f"Opening alignment file {bam}")
    alnfile = pysam.AlignmentFile(bam, "rb")
    # Check for BAM index
    try:
        alnfile.check_index()
    except:
        # check_index() should return ValueError if no index found
        logger.error(f"No valid BAM index found for alignment file {bam}")
    # Report summary to log
    logger.info("Alignment file contains "
                + str(alnfile.mapped)
                + " reads mapped to "
                + str(alnfile.nreferences)
                + " reference sequences")

    # If Fasta file provided
    if fasta:
        # Read reference Fasta file into memory
        logger.info(f"Reading mapping sequence file to memory {fasta}")
        refgenome = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

        # Sanity check: lengths and names of contigs in refgenome should match BAM
        for ctg in refgenome:
            if ctg in alnfile.references:
                if len(refgenome[ctg]) != alnfile.get_reference_length(ctg):
                    raise Exception(
                        f"Contig {ctg} length does not match BAM header")
            else:
                logger.warn(
                    f"Contig {ctg} in reference Fasta but not in BAM header")

        # Return alignment and reference objects
        return(alnfile, refgenome)

    else:
        # If no Fasta file, return alignment object only
        return(alnfile)


def milraa(args):
    logger = logging.getLogger("main.milraa")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES MILRAA")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))
    # Read input files
    alnfile, refgenome = read_bam_ref(args.bam, args.ref)
    # Initialize new IesRecords object to store putative IESs
    iesrecords = Milraa.IesRecords(alnfile, "bam", refgenome)
    # Process alignment to find putative IESs
    logger.info("Processing alignment to find putative IESs")
    if args.contig:
        logger.info(f"Extracting only reads mapping to contig {args.contig}")
        if args.start or args.stop:
            logger.info(f"from {str(args.start)} : {str(args.stop)}")
    iesrecords.findPutativeIes(args.min_ies_length,
                               args.contig, args.start, args.stop)
    if args.dump:
        logger.info(f"Dumping data in JSON format to file {args.out}.dump")
        with open(f"{args.out}.dump", "w") as fh:
            # sys.stderr.write(str(iesrecords) + "\n") # Print summary of IesRecords object
            fh.write(iesrecords.dump())  # Dump data to check
    # Get consensus locations and IES sequences
    if args.type == "ccs":
        logger.info(
            "Finding coordinates and consensus sequences of IESs from CCS read mappings")
        (iesgff, iesseq) = iesrecords.reportPutativeIes(args.min_break_coverage,
                                                        args.min_del_coverage)
    elif args.type == "subreads":
        logger.info(
            "Finding coordinates and consensus sequences of IESs from subread mappings")
        (iesgff, iesseq, iesgff_prob, iesseq_prob) = iesrecords.reportPutativeIesInsertSubreads(
            args.min_break_coverage, args.min_del_coverage,
            args.subreads_flank_len, args.subreads_pos_max_cluster_dist,
            args.subreads_cons_len_threshold)

    # Output files
    # Write gff version header and command line as comment
    with open(f"{args.out}.milraa_ies.gff3", "w") as fh:
        logger.info(
            f"Reporting putative IESs in GFF format to file {args.out}.milraa_ies.gff3")
        fh.write("##gff-version 3\n")
        fh.write("# " + " ".join(sys.argv) + "\n")
        # Write each GFF entry as a tab-separated line
        iesgff.gff2fh(fh, header=False)
    # Write Fasta file of putative IES sequences
    logger.info(
        f"Reporting consensus sequences of putative IESs to Fasta file {args.out}.milraa_ies.fasta")
    SeqIO.write(iesseq.values(), f"{args.out}.milraa_ies.fasta", "fasta")
    # Additional output files for subreads
    if args.type == "subreads":
        logger.info(
            f"Reporting putative problematic IESs in GFF format to file {args.out}.milraa_ies.problem.gff3")
        with open(f"{args.out}.milraa_ies.problem.gff3", "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write("# " + " ".join(sys.argv) + "\n")
            # Write each GFF entry as a tab-separated line
            iesgff_prob.gff2fh(fh, header=False)
        logger.info(
            f"Reporting consensus sequences of putative problematic IESs to Fasta file {args.out}.milraa_ies.fasta")
        SeqIO.write(iesseq_prob.values(),
                    f"{args.out}.milraa_ies.problem.fasta", "fasta")

    # Report junction sequences
    if args.junction_flank:
        junctionseqs = Milraa.getIndelJunctionSeqs(iesgff, iesseq,
                                                   refgenome, args.junction_flank)
        logger.info(
            f"Reporting flanking sequences of putative IESs to file {args.out}.junction.out")
        with open(f"{args.out}.junction.out", "w") as fh:
            fh.write("\t".join(["id",
                                "contig", "start", "end",
                                "leftflank", "rightflank", "pointer",
                                "tastart", "taend", "tapointer",
                                "ppstart", "ppend", "pppointer",
                                "indel", "ref"]) + "\n")  # header
            for junc in junctionseqs:
                fh.write("\t".join(junc) + "\n")

    if args.fuzzy_ies and args.type == "ccs":
        logger.info(
            "Finding coordinates and consensus sequences of IESs from CCS read mappings with allowance for unequal insert lengths")
        (fuzzygff, fuzzyiesseq) = iesrecords.reportPutativeIesInsertFuzzy(
            args.min_break_coverage,
            args.min_del_coverage,
            args.cluster_dist)
        logger.info(
            f"Reporting putative IESs with fuzzy lengths in GFF format to file {args.out}.milraa_ies_fuzzy.gff3")
        with open(f"{args.out}.milraa_ies_fuzzy.gff3", "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write("# " + " ".join(sys.argv) + "\n")
            fuzzygff.gff2fh(fh, header=False)
        logger.info(
            f"Reporting consensus sequences of putative fuzzy IESs to Fasta file {args.out}.milraa_ies_fuzzy.fasta")
        SeqIO.write(fuzzyiesseq.values(),
                    f"{args.out}.milraa_ies_fuzzy.fasta",
                    "fasta")

    # Close AlignmentFile
    alnfile.close()
    logger.info("Finished MILRAA")


def miser(args):
    logger = logging.getLogger("main.miser")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES MISER")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))

    # Read input files
    alnfile, refgenome = read_bam_ref(args.bam, args.ref)
    # Initialize new IesRecords object to store putative IESs
    iesrecords = Milraa.IesRecords(alnfile, "bam", refgenome)

    # Read in IES GFF file produced by Milraa
    logger.info("Reading GFF file containing putative IESs")
    iesgff = SharedFunctions.Gff()
    iesgff.file2gff(args.gff)

    # Compare mean mismatch percentage of reads with and without putative IES
    # for each putative IES
    logger.info("""
    Reporting possibly spurious IESs due to misassembly or mapped paralogs
    """)
    out_gff_split = defaultdict(
        list)  # dict to hold split GFF file keyed by diagnosis
    args.out.write("\t".join(['ID',
                              'mean_mismatch_pc_with_indel',
                              'mean_mismatch_pc_no_indel',
                              'stdev_with_indel',
                              'stdev_no_indel',
                              'statistic',
                              'p-value',
                              'num_reads_with_indel',
                              'num_reads_no_indel',
                              'diagnosis']))
    args.out.write("\n")
    for bpid in iesgff:
        ins_mm, non_mm = iesrecords.reportIndelReadMismatchPc(
            iesgff.getValue(bpid, 'seqid'),
            int(iesgff.getValue(bpid, 'start')),
            int(iesgff.getValue(bpid, 'end')),
            int(iesgff.getAttr(bpid, 'IES_length')),
            int(args.min_ies_length)
        )
        # Perform test of mismatch % if more than 2 reads with inserts
        # (otherwise stdev meaningless)
        if len(ins_mm) > 2 and len(non_mm) > 2:
            if args.spurious_ies_test == 'mann-whitney':
                # Mann-Whitney U test for whether mismatch % with indel of interest
                # is greater than without
                mwstat, mwpval = mannwhitneyu(
                    ins_mm, non_mm, alternative='greater')
            else:
                # Ward's t-test (non-equal population variances)
                mwstat, mwpval = ttest_ind(ins_mm, non_mm, equal_var=False)
            # Report
            outarr = [bpid,
                      round(stats.mean(ins_mm), 2),
                      round(stats.mean(non_mm), 2),
                      round(stats.stdev(ins_mm), 2),
                      round(stats.stdev(non_mm), 2),
                      round(mwstat, 2),
                      '%.2E' % mwpval,  # scientific notation
                      len(ins_mm),
                      len(non_mm)]
            # Diagnosis
            diagnosis = "ok"
            if stats.mean(ins_mm) > 5 or stats.mean(non_mm) > 5:
                diagnosis = "high_error"
            if len(ins_mm) > len(non_mm):
                diagnosis = 'misassembly'
            # PVAL_UNCORR = 0.05 # TODO magic number
            pval_corr = args.spurious_ies_pvalue / \
                len(iesgff)  # Bonferroni correction
            if mwpval < pval_corr and stats.mean(ins_mm) > stats.mean(non_mm):
                diagnosis = "paralog"
        elif len(non_mm) < 1:
            outarr = [bpid,
                      round(stats.mean(ins_mm), 2),
                      "NA",
                      "NA",
                      "NA",
                      "NA",
                      "NA",
                      len(ins_mm),
                      len(non_mm)]
            diagnosis = "misassembly"  # or scrambling
        else:
            outarr = [bpid,
                      round(stats.mean(ins_mm), 2),
                      round(stats.mean(non_mm), 2),
                      "NA",
                      "NA",
                      "NA",
                      "NA",
                      len(ins_mm),
                      len(non_mm)]
            diagnosis = "low_coverage"
        outarr.append(diagnosis)
        # Split the input GFF entries into each diagnosis group
        out_gff_split[diagnosis].append(iesgff.getEntry(bpid))
        # Write output
        args.out.write("\t".join([str(i) for i in outarr]))
        args.out.write("\n")
        # args.out.write(" ".join([str(i) for i in ins_mm]) + "\n")
        # args.out.write(" ".join([str(i) for i in non_mm]) + "\n")

    # Output split GFF files
    if args.split_gff:
        logger.info("Splitting input GFF entries into inferred categories")
        for diag in out_gff_split:
            outfile = f"{args.gff}.{diag}.gff3"
            with open(outfile, "w") as fh_spl:
                # Write gff version header and comment some info on this file
                fh_spl.write("##gff-version 3\n")
                fh_spl.write("# " + " ".join(sys.argv) + "\n")
                fh_spl.write(
                    f"# BleTIES MISER putative IESs classified as {diag}\n")
                for line in out_gff_split[diag]:
                    fh_spl.write("\t".join([str(i) for i in line]))
                    fh_spl.write("\n")

    # Close alignment filehandle
    alnfile.close()
    logger.info("Finished MISER")


def milret(args):
    logger = logging.getLogger("main.milret")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES MILRET")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))

    # Read BAM file - SAM not supported because we need random access
    alnfile = read_bam_ref(args.bam)

    # Initialize IesRetentionsMacOnly object
    iesretentions = Milret.IesRetentionsMacOnly(args.ies, alnfile)

    # Count mapping operations per site
    logger.info(
        f"Counting IES+ and IES- forms at each junction in file {args.ies}")
    iesretentions.findMappingOps()

    # Report retention scores to file
    logger.info("Calculating retention scores per junction")
    if args.use_ies_lengths:
        logger.info(
            f"Counting only inserts matching defined IES lengths to threshold +/- {str(args.length_threshold)}")
        iesretentions.calculateRetentionScoresMatchLengths(
            args.length_threshold)
    else:
        logger.info(
            "Counting all inserts at junctions as potential IESs, regardless of length")
        iesretentions.calculateRetentionScores()
    iesretentions.reportRetentionScores(args.out)
    if args.dump:
        iesretentions.dump(f"{args.out}.dump.json")

    alnfile.close()
    logger.info("Finished MILRET")


def milcor(args):
    logger = logging.getLogger("main.milcor")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES MILCOR")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))

    # Read BAM file - SAM not supported because we need random access
    alnfile = read_bam_ref(args.bam)

    logger.info(
        f"Counting per-read presence of IESs defined in file {args.ies}")
    iescorr = Milcor.IesCorrelationsByRead(args.ies, alnfile)
    if args.contig:
        logger.info(f"Extracting only reads mapping to contig {args.contig}")
        if args.start or args.stop:
            logger.info(f"from {str(args.start)} : {str(args.stop)}")
    if args.use_ies_lengths:
        logger.info(
            f"Counting only inserts matching defined IES lengths to threshold +/- {str(args.length_threshold)}")
    iescorr.countIesCooccurrences(args.use_ies_lengths, threshold=args.length_threshold,
                                  fetch_ctg=args.contig, fetch_start=args.start, fetch_stop=args.stop)

    out_perread_table = iescorr.summarizePerRead()
    logger.info(f"Writing output to file {args.out}.milcor.tsv")
    with open(f"{args.out}.milcor.tsv", "w") as fh:
        fh.write("\t".join(['qname', 'rname', 'start',
                            'end', 'ies_present', 'ies_absent']))
        fh.write("\n")
        for line in out_perread_table:
            fh.write("\t".join([str(i) for i in line]))
            fh.write("\n")

    if args.dump:
        logger.info(
            f"Dumping internal data to file {args.out}.milcor.dump.json for troubleshooting")
        iescorr.dump(f"{args.out}.milcor.dump.json")

    if args.bin:
        logger.info(
            f"Binning reads to likely MAC and MIC reads with threshold {str(args.bin_threshold)}")
        if args.contig:
            logger.info(
                f"Extracting only reads mapping to contig {args.contig}")
            if args.start or args.stop:
                logger.info(f"from {str(args.start)} : {str(args.stop)}")
        logger.info(f"Writing output to files:")
        logger.info(f"  {args.out}.milcor_bin_MAC.fasta")
        logger.info(f"  {args.out}.milcor_bin_MIC.fasta")
        logger.info(f"  {args.out}.milcor_bin_other.fasta")
        logger.info(f"  {args.out}.milcor_bin_noies.fasta")
        iescorr.binReads(f"{args.out}.milcor_bin_MAC.fasta",
                         f"{args.out}.milcor_bin_MIC.fasta",
                         f"{args.out}.milcor_bin_other.fasta",
                         f"{args.out}.milcor_bin_noies.fasta",
                         args.bin_threshold, args.contig, args.start, args.stop)

    alnfile.close()
    logger.info("Finished MILCOR")


def miltel(args):
    logger = logging.getLogger("main.miltel")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES MILTEL")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))

    alnfile, refgenome = read_bam_ref(args.bam, args.ref)

    logger.info("Getting softclipped sequences from aligned reads")
    cbscalls = Miltel.Miltel(alnfile, refgenome)  # initialize Miltel object
    if args.contig:
        logger.info(f"Processing only reads mapping to contig {args.contig}")
        if args.start or args.stop:
            logger.info(f"from {str(args.start)} : {str(args.stop)}")
    cbscalls.get_softclips(args.contig, args.start, args.stop)

    logger.info(
        f"Searching for telomere repeats of {args.telomere} at least {str(args.min_telomere_length)} bp long, using NCRF")
    cbscalls.find_telomeres(args.telomere, args.min_telomere_length)
    cbscalls_gff = cbscalls.report_CBS_GFF()

    logger.info(
        f"Searching for telomeres in reference contigs file {args.ref}")
    cbscalls.find_telomeres_ref(
        args.contig, args.telomere, args.min_telomere_length)
    cbscalls_refgff = cbscalls.report_ref_telomeres()
    # Combine into CBS annotations
    cbscalls_gff.combineGff(cbscalls_refgff)

    logger.info(
        f"Writing chromosome breakage sites to file {args.out}.miltel.telomeric.gff3")
    cbscalls_gff.gff2file(f"{args.out}.miltel.telomeric.gff3")

    if args.other_clips:
        logger.info(
            f"Tallying up other clipping junctions with clipped sequences at least {str(args.min_clip_length)} bp")
        othercalls_gff, othercalls_seqs = cbscalls.report_other_clips_GFF_fasta(
            args.min_clip_length)

        logger.info(
            f"Writing other clipping junctions to file {args.out}.miltel.othercalls.gff3")
        othercalls_gff.gff2file(f"{args.out}.miltel.othercalls.gff3")

        logger.info(
            f"Writing consensus seqs for other clipped sequences to file {args.out}.miltel.othercalls.fasta")
        SeqIO.write(othercalls_seqs,
                    f"{args.out}.miltel.othercalls.fasta", "fasta")

    if args.dump:
        logger.info(
            f"Dumping internal data to file {args.out}.dump.json for troubleshooting")
        with open(f"{args.out}.dump.json", "w") as fh:
            fh.write(cbscalls.dump())

    alnfile.close()
    logger.info("Finished MILTEL")


def insert(args):
    logger = logging.getLogger("main.insert")
    logger.info(f"BleTIES {__version__}")
    logger.info("Started BleTIES Insert")
    logger.info("Command line:")
    logger.info(" ".join(sys.argv))

    logger.info("Reading input files")
    refgenome = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
    gff = SharedFunctions.Gff()
    gff.file2gff(args.ies)

    if re.match(r"ins", args.mode):
        ies = SeqIO.to_dict(SeqIO.parse(args.iesfasta, "fasta"))
        logger.info(
            f"Inserting IESs to MAC reference to make MAC+IES hybrid reference")
        ins = Insert.Insert(refgenome, gff, ies)
        newrefgenome, newgff = ins.reportInsertedReference()
        outfasta = f"{args.out}.iesplus.fasta"
        outgff = f"{args.out}.iesplus.gff"
    elif re.match(r"del", args.mode):
        logger.info(
            f"Removing IESs from MAC+IES reference to make MAC-IES reference")
        dels = Insert.Insert(refgenome, gff, None)
        newrefgenome, newgff = dels.reportDeletedReference()
        outfasta = f"{args.out}.iesminus.fasta"
        outgff = f"{args.out}.iesminus.gff"

    logger.info(f"Writing output files to {outfasta}, {outgff}")
    with open(outfasta, "w") as fh:
        SeqIO.write(list(newrefgenome.values()), fh, "fasta")
    newgff.gff2file(outgff, header=True)
    logger.info("Finished Insert")
