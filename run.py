'''
-----------------------------------------------
Copyright 2013 Wellcome Trust Sanger Institute
Written by Zhihao Ding (zd1@sanger.ac.uk)
Released under the GPL
-----------------------------------------------

Created on 26 May 2013

'''

import sys
import os
import re
import pdb
import telseq.core
from datetime import datetime 

import logging as LG
LG.basicConfig(level=LG.INFO)
import argparse

LOG_LEVELS = {'debug': LG.DEBUG,
              'info': LG.INFO,
              'warning': LG.WARNING,
              'error': LG.ERROR
              }

VERSION="0.0.1"

def _getbams(bamindex):
    
    ofh = file(bamindex, 'rb')
    userautoID=False
    bams = []
    allids = []
    sid = 0
    for line in ofh:
        thisline = line.strip().split()
        if len(thisline)==1 and not userautoID: # found only 1 col
            userautoID = True # switch to auto id
            LG.warn(("No ID specified for bam %s \n"
            "Use serial number generated automatically as IDs" 
            "")%(line.strip()))
        bamfile=None
        if userautoID:
            sid += 1
            allids.append(sid)
            bamfile=thisline[0]
        else:
            allids.append(thisline[0])
            bamfile=thisline[1]

        bams.append(_check_bam(bamfile))
        
    ofh.close()
    assert len(bams) != 0, "Can't get BAM file paths from the file specified"
    assert len(allids) != 0, "Can't get BAM ID from the file specified"
    assert len(bams) == len(allids), "The number of IDs and the number of BAMs do not match"
    return allids, bams
    
def _check_required_sotware(f):
    for path in os.environ['PATH'].split(os.pathsep):
        if not os.path.exists(path):
            continue
        for file in os.listdir(path):
            if os.path.basename(f) == file:
                return os.path.join(path,f)
    return False
    
def _check_bam(bamfile):
    if re.search("bam$", bamfile):
        return bamfile
    else:
        print >> sys.stderr, "It seems that the BAM files provided in the list file, \ni.e. %s \ndoes not have a .bam suffix." % (thisbam)
        sys.exit(1)
def _example_usage():
    
    return '''

example usages:
==================
    
    Analysing a single BAM
    ---------------------------------
    python run.py --name mysample1 --output-dir tests/ tests/test.1.bam
     
    Analysing a list of BAMs 
    ---------------------------------
    In this mode, paths of the BAMs should listed in a text file, i.e. 
    named bamlist in the following examples
    TelSeq can be run in the following two ways:
    (1). Analyse BAMs sequentially (may take a very long time)
    python run.py --bams tests/bamlist --output-dir tests/
     
    (2). Analyse BAMs in parallel, i.e. on a computing cluster. Each run 
    requires an index number, ranging from 1 to the total number of BAMs, 
    indicating which BAM in the bamlist to be analysed. 
    # To analyse the 1st BAM in the bamlist tests/bamlist
    python run.py --index 1  --bams tests/bamlist --output-dir tests/
    
    # And don't forget to merge them after all runs are completed. 
    python run.py --merge --bams tests/bamlist --output-dir tests/
    '''
    
def main():
    
    print "================================================"
    print "======= TelSeq %s " %VERSION
    print "================================================"
    print "\n"
    
    if not _check_required_sotware('samtools'):
        print "Required software \"samtools\" can not be found in the paths specified in the PATH environment variable"
        sys.exit(1)
    
    LEVEL = LOG_LEVELS.get("info", LG.NOTSET)
    FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    LG.basicConfig(level=LEVEL,format=FORMAT)
    
    default_outputdir = "./"
    
    usage =  "Single BAM: %(prog)s --name bamid [--output-dir out] bamfile \n"
    usage += "       Multi BAM: %(prog)s [options] \n"
    
    parser = argparse.ArgumentParser(description='TelSeq: A tool for estimating telomere length from sequence data',
                                     usage=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=_example_usage())
    parser.add_argument('--version', action='version', version="TelSeq %s"%VERSION)
    
    parser.add_argument("-l", "--log", dest="level", action="store", 
                      help="Activate logging and set logging level. Acceptable values are [debug/info/warning/error]")
    parser.add_argument("-o", "--output-dir", dest="output_directory", action="store",  default=default_outputdir,
                      help="Output directory")
    parser.add_argument("-b", "--bams", dest="bams", action="store",  default=None,
                      help="A file that contains a list of BAM files and their corresponding IDs. The first column is the ID and "
                    "the second column is BAM file path, separated by white space or tab.")
    parser.add_argument("-n", "--name", dest="bamid", action="store",   default=None,
                      help="Only required when only one BAM file is specified.")
    parser.add_argument("-i", "--index", dest="index", action="store", type=int,  default=None,
                      help="Run the i th BAM in the BAM list file specified by --bams. "
                      "This number must be an integer between 1 and the total number "
                      "of BAMs in the file.")
    parser.add_argument("-m", "--merge", action="store_true", dest="mg", default=False,
                  help="merge the results from multiple runs into a single file")
    
    parser.add_argument('bam', metavar='bamfile', type=str, nargs='?', default=None,
                   help='Path for a single BAM file.')
    
    args = parser.parse_args()
    
    if args.level != None and args.level in LOG_LEVELS.keys():
        LEVEL = LOG_LEVELS.get(args.level, LG.NOTSET)
        FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
        LG.basicConfig(level=LEVEL,format=FORMAT)
        
    
    experimental_run = True
    experimental = False
    minimum = True
    outputindresults = False
    
    print "Command: %s \n"%(' '.join(sys.argv))
    
    mode = None
    if args.bam is not None and args.bams is None:
        # this is single bam mode
        mode = 1
        if args.bamid is None:
            print "Please specified a name using --name for the BAM \n"
            sys.exit(0)
        
    if args.bam is None and args.bams is not None and args.index is None:
        # this is multi bam mode, running through all sequentially
        mode = 2
        
    if args.bam is None and args.bams is not None and args.index is not None:
        # this is multi bam mode, running only for one sample specified by the index
        mode = 3
    
    if args.bam is None and args.mg and args.bams is not None:
        # merge the resulting files
        mode = 4
        
    if mode is None:
        LG.error("Usage not recognised.\n")
        parser.print_help()
#         print _example_usage()
        sys.exit(0)
        
    outdir = args.output_directory
    if outdir== default_outputdir: 
        LG.warn("Output directory is not specified (by the --output-dir option), the current "
                "directory will be used as the output directory.")
    
    if mode == 1:
        bam= _check_bam(args.bam)
        id = args.bamid
        print ("Analysing BAM %s with name %s"%(bam,id))
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam, 
                    experimental=experimental_run)
        ts.bamparser()
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=1, id=id, experimental=experimental, minimum=minimum)
    
    if mode == 2:
        ids, bams = _getbams(args.bams)
        print "Number oF BAMs specified:%d"%(len(bams))
        for i in range(len(bams)):
            bam=bams[i]
            id=ids[i]
            ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam,
                    experimental=experimental_run)
            ts.bamparser()
            if outputindresults:
                telseq.core.SeqPattern.integrate(outdir, atComplete=False,nbams=1, id=id, experimental=experimental, minimum=minimum)
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams), experimental=experimental, minimum=minimum)
        
    if mode == 3:
        ids, bams = _getbams(args.bams)
        print "Number oF BAMs specified:%d"%(len(bams))
        assert args.index >= 1 and args.index <= len(bams), "Index specified out of range"
        idx = args.index -1
        id = ids[idx]
        bam = bams[idx]
        print ("Analysing BAM %s with name %s"%(bam,id))
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam, 
                    experimental=experimental_run)
        ts.bamparser()
        if outputindresults:
            telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=1, id=id, experimental=experimental, minimum=minimum)
        
    if mode == 4:
        ids, bams = _getbams(args.bams)
        print "Number oF BAMs to be merged:%d"%(len(bams))
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams), experimental=experimental, minimum=minimum)
            
    print "Done. "
    
if __name__ == "__main__":
    main()
    
    
    
