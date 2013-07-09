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
import glob
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
    
#     if not _check_required_sotware('samtools'):
#         print "Required software \"samtools\" can not be found in the paths specified in the PATH environment variable"
#         sys.exit(1)
    
    LEVEL = LOG_LEVELS.get("info", LG.NOTSET)
    FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    LG.basicConfig(level=LEVEL,format=FORMAT)
    
    default_outputdir = "./"
    
    usage =  "(Single BAM) %(prog)s --name bamid [--output-dir out] bamfile \n"
    usage += "       (Multi BAM) %(prog)s [options] \n"
    
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
                  help="Merge the results from multiple runs into a single file")
    parser.add_argument("-fm", "--force-merge", action="store_true", dest="fm", default=False,
                  help="Merge the results from multiple runs into a single file, including the pre-existing ones that were generated "
                  "before this batch.")
    parser.add_argument("-vm", "--force-merge-specified", action="store_true", dest="vm", default=False,
                  help="Merge the results from multiple runs into a single file for only the BAMs specified in the bamlist file.")
    
    parser.add_argument('bam', metavar='bamfile', type=str, nargs='?', default=None,
                   help='Path for a single BAM file.')
    
    group = parser.add_argument_group('Experimental options')
    group.add_argument("-e", "--experimental", action="store_true", dest="expr", default=False,
                  help="Calculate varialbles that are useful for experimental purposes. ")
    group.add_argument("-v", "--details", action="store_true", dest="det", default=False,
                  help="Only output the telomere length.")
    group.add_argument("-t", "--outputeach", action="store_true", dest="ot", default=False,
                  help="Output text table result for each run.")
    
    args = parser.parse_args()
    
    if args.level != None and args.level in LOG_LEVELS.keys():
        LEVEL = LOG_LEVELS.get(args.level, LG.NOTSET)
        FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
        LG.basicConfig(level=LEVEL,format=FORMAT)
        
    
    experimental_run = experimental = args.expr
    minimum = not args.det
    outputindresults = args.ot
    
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
    
    if args.bam is None and (args.fm or args.mg or args.vm) and args.bams is not None:
        # merge the resulting files
        mgflag= args.mg + args.fm + args.vm
        if mgflag !=1:
            mode = None
        else:
            fg = [f for f in [0,1,2] if [args.mg,args.fm,args.vm][f] == True]
            mode = 4 + fg[0]
        
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
        telseq.core.SeqPattern.integrate(outdir, ids=[id], experimental=experimental, minimum=minimum)
        
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
                telseq.core.SeqPattern.integrate(outdir, ids=[id], experimental=experimental, minimum=minimum)
        sids = _check_results(outdir,ids)
        telseq.core.SeqPattern.integrate(outdir, ids=sids, experimental=experimental, minimum=minimum)
        
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
            telseq.core.SeqPattern.integrate(outdir, ids =[id], experimental=experimental, minimum=minimum)
        
    if mode == 4:
        ids, bams = _getbams(args.bams)
        sids = _check_results(outdir,ids)
        telseq.core.SeqPattern.integrate(outdir, ids = sids, experimental=experimental, minimum=minimum)
    if mode == 5:
        ids, bams = _getbams(args.bams)
        sids = _check_results(outdir,ids, force_all=True)
        telseq.core.SeqPattern.integrate(outdir, ids = sids, experimental=experimental, minimum=minimum)
    if mode == 6:
        ids, bams = _getbams(args.bams)
        sids = _check_results(outdir,ids, force_bamlist=True)
        telseq.core.SeqPattern.integrate(outdir, ids = sids, experimental=experimental, minimum=minimum)
    
    print "Done. "
    
    
def _check_results(outdir, ids, force_all=False, force_bamlist=False):
    found = glob.glob(os.path.join(outdir,'tmp',"*.bamscan.pickle"))
    res_ids = [os.path.basename(f).split(".")[0] for f in found]
    itx = [i for i in ids if i in res_ids]
    if len(itx) == 0:
        LG.info("%d These ids are found in the bamlist \n "%(len(res_ids)))
        print ids
        LG.warn("But no matching result files found in \n %s"%outdir)
        if len(res_ids) > 0:
            LG.info("Instead, found %s"%(res_ids))
        sys.exit(1)
        
    if len(itx) == len(ids) and len(ids) == len(res_ids):
        return itx
    if force_all:
        return itx
    
    if len(itx) < len(ids):
        LG.info("BAMs specified haven't all been analysed yet.")
        LG.info("Found %d bams"%(len(itx)))
        LG.info("Specified %d bams in the bamlist file"%(len(ids)))
        LG.info("Unfinished BAMs")
        LG.info([i for i in ids if i not in itx])
        LG.info("Use flag -fm, instead of -m, to merge them together anyway.")
        sys.exit(0)
        
    if len(itx) == len(ids) and len(ids) < len(res_ids):
        if force_bamlist:
            return ids
        LG.info("Found more result files than that specified in the bamlist file.")
        LG.info("%d bams specified"%(len(ids)))
        LG.info("%d result files found"%(len(itx)))
        LG.info("Extra BAMs are")
        LG.info([i for i in res_ids if i not in ids])
        LG.info("Flag -fm can be used to merge them together anyway.")
        LG.info("Flag -vm can be used to merge only the ones that are specified in the bamlist.")
        sys.exit(0)
        
if __name__ == "__main__":
    main()
    
    
    
