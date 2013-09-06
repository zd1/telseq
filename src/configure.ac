#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Telseq], [0.0.1], [zd1@sanger.ac.uk])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Telseq/telseq.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for header files.
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])
 
# Checks for library functions.
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])

# Check for the bamtools library path
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here

#AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
#	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)])

#if test -d "$with_bamtools"; then
#    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib"
#    bamtools_include="-I$with_bamtools/include"
#fi

AC_ARG_WITH(bamtools,
        AS_HELP_STRING([--with-bamtools=DIR], [specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]),
        [ if test -d "$withval" ; then
        	bamtools_include="-I$withval/include"
			bamtools_ldflags="-Wl,-rpath,$withval/lib -L$withval/lib"
		else
          	AC_MSG_ERROR([Cannot stat directory $withval])
      	fi ] 
        )

AC_ARG_WITH(gzstream,
        AS_HELP_STRING([--with-gzstream=DIR], [specify directory containing the gzstream library (http://www.cs.unc.edu/Research/compgeom/gzstream/)]),
        [ if test -d "$withval" ; then
        	gzstream_include="-I$withval"
			gzstream_ldflags="-Wl,-rpath,$withval -L$withval"
		else
          	AC_MSG_ERROR([Cannot stat directory $withval])
      	fi ] 
        )

#with_bamtools="/Users/zd1/opt/software/bamtools-2.3.0"
#bamtools_include="-I$with_bamtools/include"
#bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib"

#with_gzstream="/Users/zd1/opt/software/gzstream"
#gzstream_include="-I$with_gzstream"
#gzstream_ldflags="-Wl,-rpath,$with_gzstream -L$with_gzstream"

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror -Wno-unknown-pragmas")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $bamtools_include $gzstream_include")
AC_SUBST(LDFLAGS, "$bamtools_ldflags $gzstream_ldflags $LDFLAGS")

AC_MSG_NOTICE(["CPPFLAGS is $CPPFLAGS"])
AC_MSG_NOTICE(["LDFLAGS is $LDFLAGS"])

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS -lbamtools")

AC_MSG_NOTICE(["LIBS is $LIBS"])

# Make sure the bamtools headers can be found
AC_CHECK_HEADER([api/BamReader.h],,[AC_MSG_ERROR(
[The bamtools library must be installed (http://github.com/pezmaster31/bamtools). \
You can specify its path with the --with-bamtools=PATH option])])

# Make sure the gzstream headers can be found
AC_CHECK_HEADER([gzstream.h],,[AC_MSG_ERROR(
[gzstream library must be installed(http://www.cs.unc.edu/Research/compgeom/gzstream/). \ 
You can specify its path with the --with-gzstream=PATH option])])

AC_CONFIG_FILES([Makefile
		Telseq/Makefile
		Util/Makefile
		Thirdparty/Makefile
		 doc/Makefile 
		 examples/Makefile 
		 man/Makefile 
		 scripts/Makefile])

AC_OUTPUT



