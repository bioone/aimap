#!/usr/bin/env python3 
#
# aimap.py
#
# This code is part of the aimap package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import json
import logging
import logging.handlers
import os
import subprocess
import pandas as pd
import gffutils
import random
import shutil
import sys
import tarfile
import time
import traceback
from argparse import ArgumentParser

from aimap import (aimap_tools, aimap_config)
from aimap import __version__ as VERSION


# Process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="aimap.py")
    parser.add_argument('--version', action='version',
                        version='%(prog)s: aimap ' + VERSION)
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None, required=True,
                        help="Output directory (required)")
    parser.add_argument("-g", "--genome", dest="genomename",
                        action="store", default=None, required=True,
                        help="Genome file name (required)")
    parser.add_argument("-a", "--annotation", dest="annotation",
                        action="store", default=None,required=True,
                        help="Genome annotation file, gff3 format")
    parser.add_argument("--outfile_name", dest="outfile_name",
                        action="store", default=None,required=True,
                        help="output file name")
    parser.add_argument( "-l","--length", dest="length",
                        action="store", default=80,
                        help="Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. Default: 80 bp.")
    parser.add_argument("-f", "--filename", dest="filename",
                        action="store", default=None,
                        help="Fastq file if model is single")
    parser.add_argument("-f1", "--filename1", dest="filename1",
                        action="store", default=None,
                        help="Fastq1 file if model is paired")
    parser.add_argument("-f2", "--filename2", dest="filename2",
                        action="store", default=None,
                        help="Fastq2 file if model is paired")
    parser.add_argument("-m", "--model", dest="model",
                        action="store", default="paired",
                        choices=["single", "paired"],
                        help="aimap model (default paired)")
    parser.add_argument("-t", "--threads", dest="threads",
                        action="store", default=1,
                        help="Number of additional threads to use (default 1)")
    parser.add_argument("--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--editing_level", dest="editing_level",
                        action="store", default=0.03, 
                        help="A-I editing_level,discard position that became lower than editing_level. Default: 0.03")
    parser.add_argument("--coverage", dest="coverage",
                        action="store", default=30, 
                        help="Discard position that became lower than coverage INT. Default: 30")    
    return parser.parse_args()

	# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))
											
# Create output directory if it doesn't exist
def make_outdir():
    """Make the output directory, if required.

    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:

    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        logger.error("Output directory %s would overwrite existing " +
                         "files (exiting)", args.outdirname)
        sys.exit(1)
    else:
        logger.info("Creating directory %s", args.outdirname)
        os.makedirs(args.outdirname)
			
# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline()
	
	# Set up logging
    logger = logging.getLogger('aimap.py: %s' %
                               time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
	
	# Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging",
                         args.logfile)
            sys.exit(1)
			
	# Have we got required args? If not, exit.
    if args.genomename is None:
        logger.error("No genome file name (exiting)")
        sys.exit(1)
    if args.annotation is None:
        logger.error("No annotation file (exiting)")
        sys.exit(1)
    if args.outfile_name is None:
        logger.error("No outfile name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)

    #Step1 Apply adapter and quality trimming to FastQ files.
    logger.info("Step1 Apply adapter and quality trimming to FastQ files.")
    if args.model == "single":
        if args.filename is None:
            logger.error("No fastq file name (exiting)")
            sys.exit(1)
        else:
            os.system(aimap_tools.construct_trim_galore_single_cmdline(args.outdirname,args.length,args.filename))
    if args.model == "paired":
        if args.filename1 is None:
            logger.error("No required fastq file1 name (exiting)")
            sys.exit(1)
        elif args.filename2 is None:
            logger.error("No required fastq file2 name (exiting)")
            sys.exit(1)
        else:
           os.system(aimap_tools.construct_trim_galore_paired_cmdline(args.outdirname,args.length,args.filename1,args.filename2))
							
    #Step2 Index genome sequences in the FASTA format.
    logger.info("Step2 Index sequences in the FASTA format.")
    os.system(aimap_tools.construct_bwa_index_cmdline(args.genomename))
	
    #Step3 Map reads to genome.
    logger.info("Step3 Map reads to genome.")
    if args.model == "single":
        old_file=os.path.splitext(os.path.split(args.filename)[1])[0]
        new_file=args.outdirname+"/"+old_file+"_trimmed.fq"
        os.system(aimap_tools.construct_bwa_mem_single_cmdline(args.genomename, new_file,args.outdirname,args.outfile_name,args.threads))
    if args.model == "paired":
        old_file1=os.path.splitext(os.path.split(args.filename1)[1])[0]
        old_file2=os.path.splitext(os.path.split(args.filename2)[1])[0]
        new_file1=args.outdirname+"/"+old_file1+"_val_1.fq"
        new_file2=args.outdirname+"/"+old_file2+"_val_2.fq"
        os.system(aimap_tools.construct_bwa_mem_paired_cmdline(args.genomename, new_file1,new_file2,args.outdirname,args.outfile_name,args.threads))    
	
    #Step4 Convert sam file to bam file.
    logger.info("Step4 Convert sam file to bam file.")
    insamfile=args.outdirname+"/"+args.outfile_name+".sam"
    os.system(aimap_tools.construct_samtools_view_cmdline(insamfile,args.outdirname,args.outfile_name,args.threads))
	
    #Step5 Sort the bam file.
    logger.info("Step5 Sort the bam file.")
    inbamfile=args.outdirname+"/"+args.outfile_name+".bam"
    os.system(aimap_tools.construct_samtools_sort_cmdline(inbamfile,args.outdirname,args.outfile_name,args.threads))
    	
    #Step6 Generate pileup for one or multiple BAM files.
    logger.info("Step6 Generate pileup for one or multiple BAM files.")
    sorted_bamfile=args.outdirname+"/"+args.outfile_name+"-sorted.bam"
    os.system(aimap_tools.construct_samtools_mpileup_cmdline(sorted_bamfile,args.genomename,args.outdirname,args.outfile_name,args.threads))
    
    #Step7 Convert pileup file to table.
    logger.info("Step7 Convert pileup file to table.")
    pileup_file=args.outdirname+"/"+args.outfile_name+".pileup"
    aimap_tools.convert_pileup_to_table(pileup_file,args.outdirname,args.outfile_name)
     
    #Step8 Get a_i position.
    logger.info("Step8 Get a_i position.")
    all_info_file=args.outdirname+"/"+args.outfile_name+"_all_position_info.txt"
    aimap_tools.get_a_i_map(all_info_file,args.outdirname,args.outfile_name,args.coverage,args.editing_level)
    	
    #Step9 Confirm whether the amino acid has changed.
    logger.info("Step9 Confirm whether the amino acid has changed.")
    a_i_file=args.outdirname+"/"+args.outfile_name+"_a_i.txt"
    aimap_tools.get_founction_Prokaryote(a_i_file,args.outdirname,args.outfile_name,args.genomename,args.annotation)

    		
    # Report that we've finished
    logger.info("Done: %s.", time.asctime())
    logger.info("Time taken: %.2fs", (time.time() - t0))
