#!/usr/bin/env python

'''
-----------------------------------------------
Copyright 2013 Wellcome Trust Sanger Institute
Written by Zhihao Ding (zd1@sanger.ac.uk)
Released under the GPL
-----------------------------------------------

Created on 26 May 2013

@author: zd1
'''

from __future__ import division, print_function

DOCLINES = __doc__.split("\n")

import os
import sys

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: GPL
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Biological
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: MacOS
"""

NAME                = 'telseq'
MAINTAINER          = "Zhihao Ding"
MAINTAINER_EMAIL    = "zd1@sanger.ac.uk"
DESCRIPTION         = DOCLINES[0]
LONG_DESCRIPTION    = "\n".join(DOCLINES[2:])
URL                 = "... "
DOWNLOAD_URL        = "... "
LICENSE             = 'GPL'
CLASSIFIERS         = [_f for _f in CLASSIFIERS.split('\n') if _f]
AUTHOR              = "Zhihao Ding et al."
AUTHOR_EMAIL        = "zd1@sanger.ac.uk"
PLATFORMS           = ["Windows", "Linux", "Mac OS-X", "Unix"]
MAJOR               = 1
MINOR               = 0
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

def setup_package():

    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)
    
    # Run build
    from distutils.core import setup
    import setuptools
    
    try:
        setup(
            name=NAME,
            maintainer=MAINTAINER,
            maintainer_email=MAINTAINER_EMAIL,
            description=DESCRIPTION,
            long_description=LONG_DESCRIPTION,
            url=URL,
            download_url=DOWNLOAD_URL,
            license=LICENSE,
            classifiers=CLASSIFIERS,
            author=AUTHOR,
            author_email=AUTHOR_EMAIL,
            platforms=PLATFORMS,
            version=VERSION, 
            packages = setuptools.find_packages(".") )
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return
    

if __name__ == '__main__':
    setup_package()
