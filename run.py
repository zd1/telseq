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


def now():
    timeformat="%Y%m%d_%H%M%S"
    return datetime.strftime(datetime.now(), timeformat)

def autoNaming(bams):
    names = []
    for b in bams:
        n = ""
        if re.search(r"\\",b):
            n = b.split("\\")[-1]
        elif re.search(r"/",b):
            n = b.split("/")[-1]
        else:
            n="telseq"
        n += "_"+now()
        names.append(n)
    return names

def main(argv):

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
    parser.add_option("-n", "--names", dest="bamids", action="store", type="str", default=None,
                      help="The names that match each of the BAM files provided in -b option. BAM file path would be used as name "
                      "if this option is not provided.")
    parser.add_option("-i", "--index", dest="index", action="store", type="int", default=None,
                      help="Run the i th BAM in the list of BAMs provided in -i. If not specified, all BAMs "
                      "will be analysed sequentially. This number must be an integer between 1 and the total number "
                      "of BAMs provided.")
    
    (options, args) = parser.parse_args(argv)
    
    if len(args)==0:
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
    multi=False
    if options.bams != None:
        ofh = file(options.bams,'rb')
        for line in ofh:
            thisbam=line.strip()
            if re.search("bam$", thisbam):
                bams.append(thisbam)
            else:
                print >> sys.stderr,"In the BAM list file provided, \n %s \n does not have a .bam suffix."%(thisbam)
                sys.exit(1)
        ofh.close()
        
        if options.bamids != None:
            ofh = file(options.bamids,'rb')
            for line in ofh:
                thisbamid=line.strip()
                ids.append(thisbamid)
            ofh.close()
        else:
            ids=autoNaming(bams)
        
        if options.index != None:
            assert options.index >=1 and options.index <= len(bams), "index provided is not valid. It has to be within range [%d, %d]"%(1,len(bams)) 
            index=options.index -1 # convert to be 0-based
        multi = True
    else:
        bams=[args[0]]
        if options.bamids != None:
            ids = [options.bamids]
        else:
            ids=autoNaming(bams)
    
    # run stuff
    outdir = options.output_directory
    if outdir== default_outputdir:
        LG.warn("Output directory is not specified (by the -o option), use the current"
                "directory as the output directory.")
        
    if not multi or  (multi and index is not None):
        if index is not None:
            bam = bams[index]
            id = ids[index]
        else:
            bam = bams[0]
            id = ids[0]
        ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam)
        ts.bamparser()
        
    else:
        for i in range(len(bams)):
            bam=bams[i]
            id=ids[i]
            ts = telseq.core.SeqPattern(
                    outdir=outdir,
                    name=id,
                    bamfile=bam)
            ts.bamparser()
    
    telseq.core.SeqPattern.integrate(outdir)

if __name__ == "__main__":
    main(argv=sys.argv[1:])
    
    
    