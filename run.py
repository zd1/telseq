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
from optparse import OptionParser

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
    assert len(bams) != 0, "Can't get bam file paths from the file specified"
    assert len(allids) != 0, "Can't get bam ID from the file specified"
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

def main(argv):
    
    print "================================================"
    print "======= Telseq v%s " %VERSION
    print "================================================"
    
    if not _check_required_sotware('samtools'):
        print "Required software \"samtools\" can not be found in the paths specified in the environment PATH variable"
        sys.exit(1)
    
    LEVEL = LOG_LEVELS.get("info", LG.NOTSET)
    FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    LG.basicConfig(level=LEVEL,format=FORMAT)
    
    default_outputdir = "./"
    
    usage =  "Usage: %prog [options] [-n name] bamfile \n"
    usage += "Usage: %prog [options] [-i index] [-n bamIDs] -b bamlist "
    
    parser = OptionParser(usage=usage)
    
    parser.add_option("-l", "--log", dest="level", action="store", type="string",
                      help="Activate logging and set logging level. Acceptable values are [debug/info/warning/error]")
    parser.add_option("-o", "--output-dir", dest="output_directory", action="store", type="str", default=default_outputdir,
                      help="Output directory")
    parser.add_option("-b", "--bams", dest="bams", action="store", type="str", default=None,
                      help="A file that contains a list of BAM files, each in a separate row.")
    parser.add_option("-n", "--name", dest="bamid", action="store", type="str", default=None,
                      help="Only required when only one BAM file is specified.")
    parser.add_option("-i", "--index", dest="index", action="store", type="int", default=None,
                      help="Run the i th BAM in the list of BAMs provided in -i. If not specified, all BAMs "
                      "will be analysed sequentially. This number must be an integer between 1 and the total number "
                      "of BAMs provided.")
    parser.add_option("-m", "--merge", action="store_true", dest="mg", default=False,
                  help="merge data for multiple bams into a single file")
    
    (options, args) = parser.parse_args(argv)
    
    if options.level != None and options.level in LOG_LEVELS.keys():
        LEVEL = LOG_LEVELS.get(options.level, LG.NOTSET)
        FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
        LG.basicConfig(level=LEVEL,format=FORMAT)
        
    outdir = options.output_directory
    if outdir== default_outputdir: 
        LG.warn("Output directory is not specified (by the -o option), use the current"
                "directory as the output directory.")
    
    experimental = True
    outputindresults = True
    
    
    print "\tCommand:run.py %s"%(" ".join(argv))

    mode = None
    if len(args)==1 and options.bams is None:
        # this is single bam mode
        mode = 1
        if options.bamid is None:
            print "Please specified a name using --name for the BAM \n"
            sys.exit(0)
        
    if len(args)==0 and options.bams is not None and options.index is None:
        # this is multi bam mode, running through all sequentially
        mode = 2
        
    if len(args)==0 and options.bams is not None and options.index is not None:
        # this is multi bam mode, running only for one sample specified by the index
        mode = 3
    
    if len(args)==0 and options.mg and options.bams is not None:
        # merge the resulting files
        mode = 4
        
    if mode is None:
        LG.warn("Usage not recognised.")
        parser.print_help()
        sys.exit(0)
    
    
    if mode == 1:
        bam= _check_bam(args[0])
        id = options.bamid
        print ("\tAnalysing BAM %s with name %s"%(bam,id))
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam, 
                    experimental=experimental)
        ts.bamparser()
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=1, id=id, experimental=experimental)
    
    if mode == 2:
        ids, bams = _getbams(options.bams)
        print "\tNumber oF BAMs specified:%d"%(len(bams))
        for i in range(len(bams)):
            bam=bams[i]
            id=ids[i]
            ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam,
                    experimental=experimental)
            ts.bamparser()
            if outputindresults:
                telseq.core.SeqPattern.integrate(outdir, atComplete=False,nbams=1, id=id, experimental=experimental)
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams), experimental=experimental)
        
    if mode == 3:
        ids, bams = _getbams(options.bams)
        print "\tNumber oF BAMs specified:%d"%(len(bams))
        assert options.index >= 1 and options.index <= len(bams), "Index specified out of range"
        idx = options.index -1
        id = ids[idx]
        bam = bams[idx]
        print ("\tAnalysing BAM %s with name %s"%(bam,id))
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam, 
                    experimental=experimental)
        ts.bamparser()
        if outputindresults:
            telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=1, id=id, experimental=experimental)
        
    if mode == 4:
        ids, bams = _getbams(options.bams)
        print "\tNumber oF BAMs to be merged:%d"%(len(bams))
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams), experimental=experimental)
            
    
if __name__ == "__main__":
    main(argv=sys.argv[1:])
    
    
    
