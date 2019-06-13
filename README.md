# qc3C - Quality control for Hi-C DNA sequencing libraries

qc3C attempts to provide a means of assessing the proportion of "signal" within a Hi-C sequencing library. That is, read-pairs which are a true product of proximity ligation, rather than self-religation or shotgun noise. To accomplish this, two modes of analysis are available: BAM analysis and assembly-free kmer analysis. 

###Installation

Installing qc3C can be accomplished directly from Github using Pip or installation can be avoided by use of Docker or Singularity. One complication of the Pip based installation is that the supporting tool Jellyfish, which also provides a required Python module, must be installed separately.

####Installing Jellyfish

From a working directory of your choosing:

1. pull down the latest source from Github
2. extract and enter the project directory
3. extend the pkgconf path variable to make sure the required packing information is found at build time.
4. build and install the command line `jellyfish` tool, which is requied to create k-mer libraries.
5. enter the python hook subdirectory
6. explicitly build and install the Python module (dna_jellyfish).

All of these steps should complete without error.

_Note: Since qc3C is written exclusively for Python 3, we explicitly invoke a py3 interpreter `python3`. Be mindful of whether, on your system, using the simpler invocation `python` unintentionally invokes a py2 interpreter._

```bash
# Step 1
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz 

# Step 2
tar xzf jellyfish-2.2.10.tar.gz && cd jellyfish-2.2.10

# Step 3
export PKG_CONFIG_PATH=$PWD:$PKG_CONFIG_PATH

# Step 4
./configure && make && make install

# Step 5
cd swig/python

# Step 6
python3 setup.py build && python3 setup.py install
```

####Installing qc3C using Pip
```bash
pip3 install git+https://github.com/cerebis/qc3C
```


####Using Docker the image

We maintain an update to date image of qc3C on Dockerhub. The image also includes a working installation of Jellyfish. Both a minimal Alpine-based image and a Centos6-based image for older systems are available.

First, pull the required image from Dockerhub.
```bash
# most systems should first try
docker pull cerebis/qc3c:alpine

# if when you run qc3C, you receive a "kernel too old" error try instead
docker pull cerebis/qc3c:centos6
```

After successfully obtaining the image, qc3C or Jellyfish can be run as follows (_using the alpine image here_) :
```bash
# show qc3C help
docker run cerebis/qc3c:alpine qc3C -h

#  show jellyfish --help
docker run cerebis/qc3c:alpine jellyfish
```

###Using qc3C

####Bam mode
The bam analysis mode requires that Hi-C reads are first mapped to a reference, preferably the same genome as was used in producing the Hi-C library. The reference can be in the form of a closed genome or assembly contigs. We recommend using [BWA MEM](https://github.com/lh3/bwa) for this purpose, with options `-5SP`.

To produce a query-name sorted bam of Hi-C reads mapped to a chosen reference

```$bash
bwa mem -5SP contigs.fa hic_reads.fq | samtools view -F 0x904 -bS - | samtools sort -n -o hic_to_ref.bam -
```

####K-mer mode 
As a reference sequence is not necessary available at QC time, a second k-mer based approach is provided. At present, qc3C relies on [Jellyfish](https://github.com/gmarcais/jellyfish) to externally generate the k-mer library and Jellyfish's Python hooks internally during analysis.

From a Hi-C read-set, a Jellyfish k-mer library of size 24 could be generated as follows:

```$bash
jellyfish count -m 24 -s 2G -C -o 24mers.jf hic_reads.fq 
```

In either case, qc3C searches for evidence of proximity ligations to infer whether or not a generated library contains strong signal. Used after performing a small assessment sequencing run, this information allows researchers to choose to a full sequencing run that is more or less deep, or abandon a library entirely. 

Example usage for a library which used both Sau3AI and MluCI restriction enzymes

**BAM mode**
```$bash
> qc3C bam --mean-insert 500 -enzyme Sau3AI --enzyme MluCI --bam hic_to_ref.bam

```
**Kmer mode**
```$bash
> qc3C kmer --mean-insert 500 -enzyme Sau3AI --reads hic_reads.fq --lib 24mers.jf
```


**Command help**

```$bash
usage: qc3C [-h] [-V] {bam,kmer} ...

qc3C: Hi-C quality control

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  Version

commands:
  Valid commands

  {bam,kmer}     choose an analysis stage for further options
```

```$bash
usage: qc3C bam [-h] [-v] [-p SAMPLE_RATE] [-s SEED] [-t N] -e NEB_NAME -m
                MEAN_INSERT [-q MIN_MAPQ] -b BAM

Alignment-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose output
  -p SAMPLE_RATE, --sample-rate SAMPLE_RATE
                        Sample only a proportion of all read-pairs [None]
  -s SEED, --seed SEED  Random seed used in sampling the read-set
  -t N, --threads N     Number of threads
  -e NEB_NAME, --enzyme NEB_NAME
                        One or more case-sensitive NEB enzyme names (Use
                        multiple times for multiple files enzymes)
  -m MEAN_INSERT, --mean-insert MEAN_INSERT
                        Mean fragment length to use in estimating the
                        unobserved junction rate
  -q MIN_MAPQ, --min-mapq MIN_MAPQ
                        Minimum acceptable mapping quality [60]
  -b BAM, --bam BAM     Input name-sorted bam file of Hi-C reads mapped to
                        references
```

```$bash
usage: qc3C kmer [-h] [-v] [-p SAMPLE_RATE] [-s SEED] [-t N] -e NEB_NAME -m
                 MEAN_INSERT [--output-table OUTPUT_TABLE] [-x MAX_COVERAGE]
                 -l KMER_LIB -r FASTQ_FILE

Kmer-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose output
  -p SAMPLE_RATE, --sample-rate SAMPLE_RATE
                        Sample only a proportion of all read-pairs [None]
  -s SEED, --seed SEED  Random seed used in sampling the read-set
  -t N, --threads N     Number of threads
  -e NEB_NAME, --enzyme NEB_NAME
                        One or more case-sensitive NEB enzyme names (Use
                        multiple times for multiple files enzymes)
  -m MEAN_INSERT, --mean-insert MEAN_INSERT
                        Mean fragment length to use in estimating the
                        unobserved junction rate
  --output-table OUTPUT_TABLE
                        Save the collected per-read statistics table to a file
  -x MAX_COVERAGE, --max-coverage MAX_COVERAGE
                        Ignore regions with more than this coverage [500]
  -l KMER_LIB, --lib KMER_LIB
                        Jellyfish kmer database
  -r FASTQ_FILE, --reads FASTQ_FILE
                        FastQ format reads used in making the kmer database
                        (Use multiple times for multiple files)
```