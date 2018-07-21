AC_DEFUN([AX_LIB_HTSLIB],
#
# Handle user hints
#
[AC_MSG_CHECKING([if HTS is wanted])
AC_ARG_WITH([htslib],
  AS_HELP_STRING([--with-htslib],
                 [search for hts in DIR/include and DIR/lib]),
 [if test "$withval" != no ; then
   AC_MSG_RESULT([yes])
   if test -d "$withval" ; then
     HTS_HOME="$withval"
   else
     HTS_HOME="/usr/local"
     AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
   fi
 else
   AC_MSG_RESULT([no])
 fi],
 [AC_MSG_RESULT([yes])
  HTS_HOME="/usr/local"
 ])


#
# Locate htlibs, if wanted
#
if test -n "${HTS_HOME}" ; then

	HTS_CPPFLAGS="-I${HTS_HOME}/include/htslib/"
	HTS_LDFLAGS="-L${HTS_HOME}/lib/ -Wl,-rpath,${HTS_HOME}/lib/"
        HTS_LIBS="-lhts -lm -lpthread"

        HTS_OLD_LDFLAGS=$LDFLAGS
        HTS_OLD_CPPFLAGS=$LDFLAGS

        LDFLAGS="$LDFLAGS ${HTS_LDFLAGS} ${HTS_LIBS}"
        CPPFLAGS="$CPPFLAGS ${HTS_CPPFLAGS}"

       # AC_LANG_SAVE
       # AC_LANG_C
	#AC_CHECK_LIB([z], [zlib])
        AC_CHECK_HEADER([sam.h], [ac_cv_sam_h=yes], [ac_cv_sam_h=no])
        AC_CHECK_LIB([hts], [hts_open], [ac_cv_libhts=yes], [ac_cv_libhts=no])
        # AC_LANG_RESTORE

        LDFLAGS="$HTS_OLD_LDFLAGS"
        CPPFLAGS="$HTS_OLD_CPPFLAGS"

        if test "$ac_cv_libhts" = "yes" && test "$ac_cv_sam_h" = "yes" ; then
                #
                # If both library and header were found, use them
                #
                AC_MSG_CHECKING([htslib])
                AC_MSG_RESULT([ok])
		AC_SUBST(HTS_CPPFLAGS)
        	AC_SUBST(HTS_LDFLAGS)
		AC_SUBST(HTS_LIBS)
                with_hts=yes
        else
                #
                # If either header or library was not found, revert and bomb
                #

                AC_MSG_CHECKING([htslib])
                AC_MSG_RESULT([failed])
                AC_MSG_ERROR([specify a valid hts installation with --with-htslib=DIR ])
        fi

fi
])
