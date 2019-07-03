# Overview
Ryuto is a tool for exact and fast transcript assembly and quantification, using network flows and a novel extension of splice-graphs.

Run
```
ryuto --help
```
to see options.

# Testdata

You may download the alignments of real data used for publication [here](http://silo.bioinf.uni-leipzig.de/thomas/ryuto_real_alignments2.tar.gz) (Size: 82GB).

# Pre-Compiled
Operation System | Version | Binary
 ---------------- | ------- | ------
Ubuntu           | 26      | [binary](https://github.com/studla/RYUTO/releases/download/1.3m-Ubuntu-26/ryuto)
Ubuntu           | 28      | [binary](https://github.com/studla/RYUTO/releases/download/1.3m-Ubuntu-28/ryuto)
Fedora           | 27      | [binary](https://github.com/studla/RYUTO/releases/download/1.3m-Fedora-27/ryuto)
Fedora           | 28      | [binary](https://github.com/studla/RYUTO/releases/download/1.3m-Fedora-28/ryuto)
Generic Linux    | -       | [binary](https://github.com/studla/RYUTO/releases/download/1.3m-Generic/ryuto)

Fedora users may have to install g++ in order for the pre-compiled binaries to work.
```
sudo dnf install gcc-c++
```

The generic version provides libstdc++ and libgcc static linked for linux distributions that cannot provide libraries supporting C11.

# Installation
Download the newest source code from: https://github.com/studla/RYUTO

The following additional libraries need to installed in order to run Ryuto:
zlib
(openmp)
boost
htslib
clp
lemon

Compile Ryuto with:
```
./configure [options]
make
```
The compiled binary can then be found in the `src` subfolder.

## ZLIB

htslib and thus Ryuto relies on zlib. If zlib is not installed on your system,
you have to install it first. Download and install zlib from (https://zlib.net/).

Use the following commands to install zlib:
```
./configure
make
make install
```

Fedora users may alternatively run:
```
sudo dnf install zlib-devel
```

Ubuntu users may alternatively run:
```
sudo apt-get install zlib1g-dev
```

You may specify a non-standard zlib installation with --with-zlib=/path/to/your/zlib to the Ryuto ./configure call.

## OpenMP

OpenMP is used for parallelization and needs to be installed in order for this feature to work.

## Install Boost

Download and install boost from (http://www.boost.org).

On Fedora you may run:
```
sudo dnf install boost-devel
```

On Ubuntu you may run:
```
sudo apt-get install libboost-all-dev
```

If installed to a non-standard path, add with --with-boost=/path/to/your/boost to the Ryuto ./configure call.

## Install htslib

Download and install htslib from (https://github.com/samtools/htslib).
Install with:
```
./configure
make
make install
```

You may alternatively use the following call if bz2 is not installed. 
```
./configure --disable-bz2 --disable-lzma
```

If installed to a non-standard path, add with --with-htslib=/path/to/your/htslib to the Ryuto ./configure call.

## Install Clp

Download and install clp. Please use the version provided in the `libraries_to_install` subfolder.
You may find details under (https://projects.coin-or.org/Clp).

Use the following call to install Clp to a custom location
```
./configure --disable-bzlib --disable-zlib --prefix=/path/to/your/clp
make
make install
```

If installed to a non-standard path, add with --with-clp=/path/to/your/clp to the Ryuto ./configure call..

## Install Lemon

Download and install lemon after clp. Please use the version provided in the `libraries_to_install` subfolder.
The provided version contains a simplified makefile. If you want to use your own version of lemon,
please make sure it is linked appropriately to CBC, CLP and COIN-UTILS.
You may find details under (http://lemon.cs.elte.hu/trac/lemon).

Install using
```
mkdir build
cd build
cmake -DLEMON_DEFAULT_LP=CLP -DCOIN_ROOT_DIR=/path/to/your/Clp ..
make
make install
```

If installed to a non-standard path, add with --with-lemon=/path/to/your/lemon to the Ryuto ./configure call.








