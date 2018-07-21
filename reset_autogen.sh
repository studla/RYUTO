#! /bin/sh


FILES=$(find . -type f -name "*.cc" | cut -c 7- | tr '\n' ' ')
sed -i "s|^\s*ryuto_SOURCES.*$|ryuto_SOURCES=$FILES|" src/Makefile.am

aclocal
autoconf
automake --add-missing --foreign


