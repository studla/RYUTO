AC_INIT([ttrans], [1.0])
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CXX

CXXFLAGS=" -O2"

AC_LANG([C++])

#zlib for HTSLIB
AX_CHECK_ZLIB

#SAM-TOOL to read Files
AX_LIB_HTSLIB

#openMP
AC_OPENMP
 
#Boost-Framework, for options and graphs
AX_BOOST_BASE([1.3],,[
	echo "Only tested for Bosst 1.3 and above!"
        exit -1
])
AX_BOOST_PROGRAM_OPTIONS
BOOST_LDFLAGS+=" -Wl,-rpath,`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`"

AX_LIB_LEMON

AX_LIB_CLP

AX_CHECK_COMPILE_FLAG([-std=c++0x], CXXFLAGS="$CXXFLAGS -std=c++0x")

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT


