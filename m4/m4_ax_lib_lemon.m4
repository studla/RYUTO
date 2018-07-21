AC_DEFUN([AX_LIB_LEMON],
#
# Handle user hints
#
[AC_MSG_CHECKING([if lemon is wanted])
AC_ARG_WITH([lemon],
  AS_HELP_STRING([--with-lemon],
                 [search for lemon in DIR/include and DIR/lib]),
 [if test "$withval" != no ; then
   AC_MSG_RESULT([yes])
   if test -d "$withval" ; then
     LEMON_HOME="$withval"
   else
     LEMON_HOME="/usr/local"
     AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
   fi
 else
   AC_MSG_RESULT([no])
 fi],
 [AC_MSG_RESULT([yes])
  LEMON_HOME="/usr/local"
 ])


#
# Locate htlibs, if wanted
#
if test -n "${LEMON_HOME}" ; then

	LEMON_CPPFLAGS="-I${LEMON_HOME}/include/"
	LEMON_LDFLAGS="-L${LEMON_HOME}/lib/ -Wl,-rpath,${LEMON_HOME}/lib/"
        LEMON_LIBS="-lemon"

        LEMON_OLD_LDFLAGS=$LDFLAGS
        LEMON_OLD_CPPFLAGS=$LDFLAGS

        LDFLAGS="$LDFLAGS ${LEMON_LDFLAGS} ${LEMON_LIBS}"
        CPPFLAGS="$CPPFLAGS ${LEMON_CPPFLAGS}"

       # TODO: really check
       # AC_LANG_SAVE
       # AC_LANG_C
       # AC_CHECK_HEADER([lemon/core.h], [ac_cv_lemon_h=yes], [ac_cv_lemon_h=no])
       # AC_CHECK_LIB([emon], [lemonReaderNextRecord], [ac_cv_liblemon=yes], [ac_cv_liblemon=no])
        # AC_LANG_RESTORE

        LDFLAGS="$LEMON_OLD_LDFLAGS"
        CPPFLAGS="$LEMON_OLD_CPPFLAGS"

       # if test "$ac_cv_liblemon" = "yes" && test "$ac_cv_lemon_h" = "yes" ; then
                #
                # If both library and header were found, use them
                #
                AC_MSG_CHECKING([lemon])
                AC_MSG_RESULT([ok])
		AC_SUBST(LEMON_CPPFLAGS)
        	AC_SUBST(LEMON_LDFLAGS)
		AC_SUBST(LEMON_LIBS)
                with_LEMON=yes
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
