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

def now():
    timeformat="%Y%m%d_%H%M%S"
    return datetime.strftime(datetime.now(), timeformat)

def getbams(bamindex):
    
    ofh = file(bamindex, 'rb')
    userautoID=False
    bams = []
    allids = []
    sid = 0
    for line in ofh:
        thisline = line.strip().split()
        if len(thisline)==1:
            userautoID = True
            LG.warn(("No ID specified for bam %s \n"
            "Use serial number generated automatically as IDs" 
            "")%(line.strip()))
        if userautoID:
            sid += 1
            allids.append(sid)
        else:
            allids.append(thisline[0])
            
        if re.search("bam$", thisline[1]):
            bams.append(thisline[1])
        else:
            print >> sys.stderr, "It seems that the BAM files provided in the list file, \ni.e. %s \ndoes not have a .bam suffix." % (thisbam)
            sys.exit(1)
    ofh.close()
    assert len(bams) != 0, "Can't get bam file paths from the file specified"
    assert len(allids) != 0, "Can't get bam ID from the file specified"
    assert len(bams) == len(allids), "The number of IDs and the number of BAMs do not match"
    return allids, bams

def check_required_sotware(f):
    for path in os.environ['PATH'].split(os.pathsep):
        if not os.path.exists(path):
            continue
        for file in os.listdir(path):
            if os.path.basename(f) == file:
                return os.path.join(path,f)
    return False
    
def main(argv):
    
    print "================================================"
    print "======= Telseq v%s " %VERSION
    print "================================================"
    
    if not check_required_sotware('samtools'):
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
    
    outdir = options.output_directory
    if outdir== default_outputdir:
        LG.warn("Output directory is not specified (by the -o option), use the current"
                "directory as the output directory.")
    
    # if runs have already completed, only applicable to the cases of analysing multiple BAMs, 
    # and user wants to merge all the results. 
    if options.mg == True:
        if options.bams != None:
            ids, bams = getbams(options.bams)
            telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams))
            print "\tCommand:run.py %s"%(" ".join(argv))
            print "\tMerging result files "
        
        sys.exit(0)
    
    # In the case of only a single BAM is specified, check if it's there
    if len(args)==0 and options.bams == None:
        print "\n"
        print "ERROR: no bam file has been specified \n\n"
        parser.print_help()
        sys.exit()
    
    if options.level != None and options.level in LOG_LEVELS.keys():
        LEVEL = LOG_LEVELS.get(options.level, LG.NOTSET)
        FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
        LG.basicConfig(level=LEVEL,format=FORMAT)
        
    bams = []
    ids = []
    index=None
    
    # In the case of multiple bams are specified in a text file, the bam index file, 
    # load all of them in. 
    if options.bams != None:
        if options.mg :
            telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=len(bams))    
            sys.exit(0)
        ids, bams = getbams(options.bams)
        
        if options.index != None:
            assert options.index >=1 and options.index <= len(bams), "index provided is not valid. It has to be within range [%d, %d]"%(1,len(bams)) 
            index=options.index -1 # convert to be 0-based
    else:
        bams=[args[0]]
        if options.bamid != None:
            ids = [options.bamid]
        else:
            ids= 1
    

    print "\tCommand:run.py %s"%(" ".join(argv))
    print "\tNumber oF BAMs specified:%d"%(len(bams))
    
    if options.bams == None or (options.bams != None and index != None): # we only analyse one BAM
        if index != None:
            bam = bams[index]
            id = ids[index]
            print ("\tAnalysing BAM %s with name %s"%(bam,id))
            print ("\tIt's the %d th BAM in the BAM list specified"%(index+1))
        else:
            bam = bams[0]
            id = ids[0]
            print ("\tAnalysing BAM %s with name %s"%(bam,id))
            
        if not os.path.exists(bam):
            print "\n"
            print "ERROR: Can't find BAM file\n%s\n"%bam
            sys.exit()
            
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam)
        ts.bamparser()
        telseq.core.SeqPattern.integrate(outdir, atComplete=True, nbams=1)
        
    else:
        for i in range(len(bams)):
            bam=bams[i]
            id=ids[i]
            ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam)
            ts.bamparser()




if __name__ == "__main__":
    main(argv=sys.argv[1:])
    
    
    
