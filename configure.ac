AC_INIT([ttrans], [1.0])
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CXX

CXXFLAGS="$CXXFLAGS -O2"
#CXXFLAGS="-O2"
#CXXFLAGS=" -ggdb"

AC_LANG([C++])

#zlib for HTSLIB
AX_CHECK_ZLIB

#SAM-TOOL to read Files
AX_LIB_HTSLIB

#openMP
AC_OPENMP
 
#Boost-Framework, for options and graphs
AX_BOOST_BASE([1.3],,[
	echo "Only tested for Boost 1.3 and above!"
        exit -1
])
AX_BOOST_PROGRAM_OPTIONS
BOOST_LDFLAGS+=" -Wl,-rpath,`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`"

AX_LIB_LEMON

AX_LIB_CLP

AX_CHECK_COMPILE_FLAG([-std=c++14], CXXFLAGS="$CXXFLAGS -std=c++14")

# Check whether static is enabled
AC_ARG_ENABLE([static], [AS_HELP_STRING([--enable-static], [enable static linking])])

# If statis is enabled, directly assign .a files to LIBS
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/work/ryuto_branches/RYUTO/libz.a $LIBS"]) AC_SUBST( [Z_LIBS], [""])])
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/work/ryuto_branches/RYUTO/libbz2.a $LIBS"])])
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/Libraries/lib/libhts.a $LIBS"]) AC_SUBST( [HTS_LIBS], ["-lm -lpthread"])])
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/Libraries/lib/libboost_program_options.a $LIBS"]) AC_SUBST( [BOOST_PROGRAM_OPTIONS_LIB], [""])])
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/Libraries/lib/libCoinUtils.a $LIBS"]) AC_SUBST( [CLP_LIBS], [""])])
AS_IF([test "x$enable_static" = xyes], [AC_SUBST([LIBS], ["/homes/biertruck/thomas/Libraries/lib/libClp.a $LIBS"])])

AC_ARG_WITH([staticcoin], AS_HELP_STRING([--with-staticcoin], [we use static coins clp]), [
AC_SUBST([LIBS], ["$withval/lib/libCoinUtils.a $LIBS"])
AC_SUBST([LIBS], ["$withval/lib/libClp.a $LIBS"])
], [])


#LDFLAGS="-static $LDFLAGS"


AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT


