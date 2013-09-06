#!/bin/sh

set -ex
./clean-auto-gen.sh
aclocal
autoconf
autoheader
automake -a

echo "./configure --with-bamtools=/Users/zd1/opt/software/bamtools-2.3.0 --with-gzstream=/Users/zd1/opt/software/gzstream"
echo "./configure --with-bamtools=/nfs/users/nfs_z/zd1/opt/software/bamtools-2.3.0 --with-gzstream=/nfs/users/nfs_z/zd1/opt/software/gzstream"
