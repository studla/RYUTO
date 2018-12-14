# Overview
Ryuto is a tool for exact and fast transcript assembly and quantification, using network flows and a novel extension of splice-graphs.

# Testing
The alignments used for testing can be downloaded here: [real datasets](http://silo.bioinf.uni-leipzig.de/thomas/ryuto_real_alignments.tar.gz) and [simulated datsets](http://bioinf.itmat.upenn.edu/BEERS/bp2/)

# Installation
Download the source code from: https://github.com/studla/RYUTO

The following additional libraries need to installed in order to run Ryuto:
boost
htslib
lemon
clp

Compile with:
./configure [options]
make

## 1. Install Boost

Download and install boost from http://www.boost.org.

If installed to a non-standard path, add with --with-boost for configure.

## 2. Install htslib

Download and install htslib from https://github.com/samtools/htslib

If installed to a non-standard path, add with --with-htslib for configure.

## 3. Install Clp

Download and install clp from https://projects.coin-or.org/Clp

If installed to a non-standard path, add with --with-clp for configure.
In this case you will also need to add -DLEMON_DEFAULT_LP=CLP -DCOIN_ROOT_DIR=custom_path

## 4. Install Lemon

Download and install lemon from http://lemon.cs.elte.hu/trac/lemon

If installed to a non-standard path, add with --with-lemon for configure.

## 5. ZLIB

You may specify a non-standard zlib installation with --with-zlib or disable it with --without-zlib .

## 6. OpenMP

OpenMP is used for parallelization and needs to be installed in order for this feature to work.



