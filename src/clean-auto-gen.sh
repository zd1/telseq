#!/bin/bash

if [ -h missing ]; then
    unlink missing
fi 

if [ -h install-sh ]; then 
    unlink install-sh
fi 

if [ -h compile ]; then 
    unlink compile
fi 

if [ -h depcomp ]; then 
    unlink depcomp
fi 

if [ -e config.h ]; then 
    rm config.h
fi

if [ -e aclocal.m4 ]; then 
    rm aclocal.m4
fi 

if [ -e configure ]; then 
    rm configure
fi 

if [ -d autom4te.cache ]; then 
    rm -r autom4te.cache
fi

if [ -e config.status ]; then 
    rm config.status
fi

if [ -e Makefile ]; then 
    rm Makefile
fi

if [ -e Makefile.in ]; then
    rm Makefile.in
fi

if [ -e stamp-h1 ]; then 
    rm stamp-h1
fi

if [ -e config.h.in ]; then 
    rm config.h.in
fi

if [ -e config.log ]; then 
    rm config.log
fi

if [ -e autoscan.log ]; then 
    rm autoscan.log
fi

if [ -e *~ ]; then 
    rm *~
fi

if [ -e configure.scan ]; then 
    rm configure.scan
fi

