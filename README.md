# Overview
Ryuto is a tool for exact and fast transcript assembly and quantification, using network flows and a novel extension of splice-graphs.

**Cite:**
*Gatter, Thomas, and Peter F. Stadler. "Ryūtō: network-flow based transcriptome reconstruction." BMC bioinformatics 20.1 (2019): 190.*


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

# Use

In its basic usecase, Ryuto needs to be provided only an output directory and the library type used by RNA-Seq: You must provide at least one mapping file in sam or bam format.

```
ryuto -o [output dir] -l [fr-unstranded | fr-firststrand | fr-secondstrand]  <input1.bam> <input2.bam> ...
```

We recommnd to use multiple threads via option `-t [CPU COUNT]`.
If an annotation is available for your organism, you may provide it via option `-g [GTF]`.

See 
```
ryuto --help
```
for a full list of options.

# Mapping

For best results in combination with Ryuto, especially for multi-sample assembly, we recommend Mapping reads with [STAR](https://github.com/alexdobin/STAR).

# Output

3 files will be produced:
- transcripts.gtf: The main result of the assembly. Transcripts are given in GTF format.
- transcripts.count: Counting table produced for e.g. differential transcript expression analysis. Formatted as a table as follows: `Gene-Name, Transcript-Name,	Length of the Transcript,	Read Count Sample 1 [, Read Count Sample 2 [, Read Count Sample 3 ...]]`
- transcripts.errcount: If multiple samples are assembled into a consensus, individual inputs may provide signals on how strongly they disagree with it. Higher numbers indicate a higher level of disagreement, e.g. because of undetected fold-changes. Formatted as a table as follows: `Chromosome, Strand, Start, End, Feature (Exon or Splice-Junction), Disagreement Sample 1 [, Disagreement Sample 2 [, Disagreement Sample 3 ...]]`

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


# Testdata

You may download the alignments of real data used for the first publication [here](http://silo.bioinf.uni-leipzig.de/thomas/ryuto_real_alignments2.tar.gz) (Size: 82GB).







