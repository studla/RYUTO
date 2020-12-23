# Overview
Ryuto is a tool for exact and fast transcript assembly and quantification, using network flows and a novel extension of splice-graphs.

Run
```
ryuto --help
```
to see options.

# Bioconda

Ryuto is available via [bioconda](https://bioconda.github.io/)

```
conda install ryuto
```

# Pre-Compiled
Operation System |  Binary
 ---------------- | ------
Generic Linux, dynamic linked (recommended) | [binary](https://github.com/studla/RYUTO/releases/download/1.6/ryuto)
Generic Linux, static libc | [binary](https://github.com/studla/RYUTO/releases/download/1.6-static/ryuto)

All users have to install g++ in order for the pre-compiled dynamic binaries to work.
e.g, for Fedora use:
```
sudo dnf install gcc-c++
```

The static version provides libstdc++ and libgcc static linked for linux distributions that cannot provide libraries supporting >= C11.

# Installation
Download the newest source code from: https://github.com/studla/RYUTO

The following additional libraries need to installed in order to run Ryuto:
zlib
openmp
boost
htslib

Ryuto provides the following packages with simplified installers.
clp
lemon

Compile Ryuto with:
```
./build_all.sh --prefix=[Install Path] [options]
```
The compiled binary is installed to the specified prefix.

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

You may specify a non-standard zlib installation with --with-zlib=/path/to/your/zlib to the Ryuto ./build_all.sh call.

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

If installed to a non-standard path, add with --with-boost=/path/to/your/boost to the Ryuto ./build_all.sh call.

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

If installed to a non-standard path, add with --with-htslib=/path/to/your/htslib to the Ryuto ./build_all.sh call.

## Install Clp and Lemon

Clp and Lemon are automatically built and during the ./build_all.sh call. You can find both in the "extern" subfolder.

Lemon needs to be linked to CLP if you want to provide those libraries yourself.
If installed to a non-standard path, add with --with-clp=/path/to/your/clp --with-lemon=/path/to/your/lemon and build the project manually.

Compile Ryuto with:
```
./configure [options]
make
```

# Testdata

You may download the alignments of real data used for the first publication [here](http://silo.bioinf.uni-leipzig.de/thomas/ryuto_real_alignments2.tar.gz) (Size: 82GB).







