#!/bin/sh

set -ex
./clean-auto-gen.sh
aclocal
autoconf
autoheader
automake -a

