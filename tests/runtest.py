#!/usr/bin/env python

'''
-----------------------------------------------
Copyright 2013 Wellcome Trust Sanger Institute
Written by Zhihao Ding (zd1@sanger.ac.uk)
Released under the GPL
-----------------------------------------------

Created on 26 May 2013

'''

PROJECT_MODULE = "telseq"
PROJECT_ROOT_FILES = ['telseq', 'LICENSE', 'setup.py']
SAMPLE_TEST = ""
SAMPLE_SUBMODULE = ""

# ---------------------------------------------------------------------

import sys
import os

os.chdir(sys.path[0])
sys.path.append("./../")

import telseq.core

import subprocess
from argparse import ArgumentParser

def main():
    
    ts = telseq.core.SeqPattern(
                    outdir="./",
                    name='test',
                    bamfile="./test.1.bam")
    ts.bamparser()
    telseq.core.SeqPattern.integrate("./", atComplete=True, nbams=1)
    
if __name__ == "__main__":
    main()

