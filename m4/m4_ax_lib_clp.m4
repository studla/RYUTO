AC_DEFUN([AX_LIB_CLP],
#
# Handle user hints
#
[AC_MSG_CHECKING([if clp is wanted])
AC_ARG_WITH([clp],
  AS_HELP_STRING([--with-clp],
                 [search for clp in DIR/include and DIR/lib]),
 [if test "$withval" != no ; then
   AC_MSG_RESULT([yes])
   if test -d "$withval" ; then
     CLP_HOME="$withval"
   else
     CLP_HOME="/usr/local"
     AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
   fi
 else
   AC_MSG_RESULT([no])
 fi],
 [AC_MSG_RESULT([yes])
  CLP_HOME="/usr/local"
 ])


#
# Locate htlibs, if wanted
#
if test -n "${CLP_HOME}" ; then

	CLP_CPPFLAGS="-I${CLP_HOME}/include/"
	CLP_LDFLAGS="-L${CLP_HOME}/lib/ -Wl,-rpath,${CLP_HOME}/lib/"
        CLP_LIBS="-lClp -lCoinUtils"

        CLP_OLD_LDFLAGS=$LDFLAGS
        CLP_OLD_CPPFLAGS=$CPPFLAGS

        LDFLAGS="$LDFLAGS ${CLP_LDFLAGS} ${CLP_LIBS}"
        CPPFLAGS="$CPPFLAGS ${CLP_CPPFLAGS}"

       # TODO: really check
       # AC_LANG_SAVE
       # AC_LANG_C
       # AC_CHECK_HEADER([lemon/core.h], [ac_cv_lemon_h=yes], [ac_cv_lemon_h=no])
       # AC_CHECK_LIB([emon], [lemonReaderNextRecord], [ac_cv_liblemon=yes], [ac_cv_liblemon=no])
        # AC_LANG_RESTORE

        LDFLAGS="$CLP_OLD_LDFLAGS"
        CPPFLAGS="$CLP_OLD_CPPFLAGS"

       # if test "$ac_cv_liblemon" = "yes" && test "$ac_cv_lemon_h" = "yes" ; then
                #
                # If both library and header were found, use them
                #
                AC_MSG_CHECKING([clp])
                AC_MSG_RESULT([ok])
		AC_SUBST(CLP_CPPFLAGS)
        	AC_SUBST(CLP_LDFLAGS)
		AC_SUBST(CLP_LIBS)
                with_CLP=yes
      #  else
                #
                # If either header or library was not found, revert and bomb
                #

       #         AC_MSG_CHECKING([LEMON])
       #         AC_MSG_RESULT([failed])
       #         AC_MSG_ERROR([specify a valid lemon installation with --with-lemon=DIR ])
       # fi

fi
])
