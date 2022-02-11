# qc3C - reference-free quality control for Hi-C DNA sequencing libraries

The primary aim of qc3C is to infer the proportion of "signal" within a Hi-C sequencing library and thereby provide a means of quality control. A big advantange of qc3C is that this can be done without access to a reference sequence, which until now has been a significant stopping point for projects not involving model organisms. 

By "signal" we mean the proportion of read-pairs which are true products of proximity ligation rather than pairs resulting from undesirable processes, which in downstream analyses will constitute noise in the experiment. When inferring the proportion of Hi-C pairs, no restriction is made on their separation distance. This is in contrast to reference-based methods which typically look at only the proportion of pairs that map at long-range (e.g. the number of pairs separated by more than 10 kbp). 

Although the proportion of long-range pairs is indeed an important attribute, this type of information cannot be established without the use of a reference sequence (preferably high-quality). In cases where a draft assembly is used, sequence fragmentation reduces the likelihood of observing pairs with large separation.

Still, details from reference-based assessment remain of interest and therefore qc3C can also perform this type of analysis. What we call "bam mode" analysis. Statistics obtained from this mode include such details as the number of read-through events and HiCPro style pair categorisation (e.g. dangling-end, self-circle, etc).

## Installation

### Requirements

Due to dependency issues, qc3C currently runs only on the Linux x86_64 platform. With that said, we maintain a Docker image which enables users to run qc3C on other platforms, both in Docker and Singularity.

### Using Conda

We maintain our own conda packages for both qc3C and a few supporting packages on the Anaconda channel [cerebis](https://anaconda.org/cerebis). 

**Note:** Older distributions of Linux (cira <=2016) may encounter the error ("glibc 2.23 not found") when attempting to run the resulting Conda environment. If this is you, we recommend that you either take advantage of our Docker image or follow the easy [two-step procedure below](#using-two-step-procedure).

Installation is accomplished as follows.
```
conda create -y -n qc3c -c cerebis -c conda-forge -c bioconda qc3C
``` 

### Please note the following regardng Jellyfish

Jellyfish is released as a separate packasge on our own Anaconda channel (`cerebis`), however at present there exists a name clash with another unrelated package within the large `conda-forge` channel. To prevent conda from obtaining the wrong package it is important to mention `-c cerebis` before `-c conda-forge`, as the order of channels on the commandline establishes priority. We will likely modify our package name to elininate this issue.

Further, do not forcibly install the alternative `bioconda` package kmer-jellyfish. qc3C requires support for Python language hooks within Jellfyish, however due to the bioconda requirement that any included package support both Linux and OSX, kmer-jellyfish has dropped all hook support since the more complex build fails under OSX. Consequently, if this package is installed instead of ours, qc3C will throw the error `NoModuleFoundError: No module named 'jellyfish'`. 

### Using Docker

A docker image is available from DockerHub at docker://cerebis/qc3c.

To use the Docker image, you must make your data available to the running Docker container through either a bind mount or shared volume. In the example below, the data is imagined to be stored within the directory `$PWD/mydata`. Access to this location is achieved through a volume, where the default working directory within the Docker container is `/app`. In doing it this way, the output directory will be created underneath `$PWD/mydata`.

```
# fetch the latest qc3C image
docker pull cerebis/qc3c

# kmer mode
docker run -v $PWD/mydata:/app cerebis/qc3c kmer --mean-insert 300 --enzyme DpnII --lib kmers.jf \
    --reads reads_r1.fq.gz --reads reads_r2.fq.gz --output-path output

# bam mode
docker run -v $PWD/mydata:/app cerebis/qc3c bam --enzyme DpnII --fasta ref.fna.gz --bam reads2ref.bam \
    --output-path output
```

### Using Singularity

We rely on singularity users importing our Docker image, which could be done as follows. 

**NOTE:** Singularity installations can automatically mount the host filesystem, which can result in masking the paths used by the qc3C container. It is recommended that users override this behaviour (`--contain`) and bind the paths that they require to access the files to be analyzed (`--bind [data-path]`).

If only the host path is specified in a bind argument, the path will be replicated in the container. To access files relative to your current working directory, simply bind `$PWD`.

```
singularity build qc3C.sif docker://cerebis/qc3c
singularity run --contain --bind $PWD qc3C.sif bam -e DpnII --fasta data/reference.fasta --bam data/hic2ref.bam
```

### Using Two-Step Procedure

This procedure employs a combination of Conda and Pip to install qc3C.

#### 1. Create a Conda environment with Jellyfish and Python

qc3C can be installed directly from github, however this approach requires that [Jellyfish](https://github.com/gmarcais/Jellyfish), along with its Python language hooks be installed first. Although the basic Jellyfish binaries are easily built and installed, a build which correctly installs the necessary Python hooks can be more problematic. [See Jellyfish issue #134](https://github.com/gmarcais/Jellyfish/issues/134) for more information on how to remedy this issue. 

Instead, we strongly encourage users to install our Jellyfish Conda package.

**Please Note** Do not use Bioconda's Jellyfish package, as it contains only the Jellyfish binaries and no language hooks. As a result, qc3C will fail to run with the error `NoModuleFoundError: No module named 'jellyfish'`. 

#### 2. Install qc3C using Pip directly from Github.
 
```$bash
# step 1
conda create -y -n qc3c -c cerebis -c conda-forge -c bioconda jellyfish
conda activate qc3c
# step 2
pip install git+https://github.com/cerebis/qc3C
```

### Known Problems Running qc3C

1. An error at runtime "glibc 2.23 not found"

    Older distributions of Linux may not possess a recent enough version of the glibc library. In these cases, we strongly urge users to use the [two-step procedure outlined above](#using-two-step-procedure) or the Docker image.


## Using qc3C

### Modes of analysis

qc3C is capable of quality assessing Hi-C read-sets in two ways. The traditional method, reliant on aligning reads to a reference sequence, and a new reference-free _k_-mer based approach.

### General Requirements

- Adapter trimmed Hi-C read-set in FastQ format.

We recommend [fastp](https://github.com/OpenGene/fastp) for adapter trimming.

#### Per-mode Requirements

- **Reference-free _k_-mer mode** 
  - mean insert size
  - the names of the restriction enzymes used in protocol digest 

- **Reference-dependent bam mode**
  - a (preferably high quality) reference sequence
  - the names of the restriction enzymes used in protocol digest

### Quick start

Below is a quick demonstration of how you could invoke qc3C in either _k_-mer or bam mode.

#### Creating a _k_-mer database and then analysing the library

For a Hi-C library created with the 4-cutter enzyme DpnII and a known average insert size of 300 bp.

1. first we create the _k_-mer database with the subcommand `mkdb` using default settings.
2. the _k_-mer mode analysis is then subsequently run
```
qc3C mkdb --reads reads_r1.fq.gz --reads reads_r2.fq.gz --lib kmers.jf
qc3C kmer --mean-insert 300 --enzyme DpnII --lib kmers.jf --reads reads_r1.fq.gz --reads reads_r2.fq.gz --output-path output
```

#### Analysing in _k_-mer mode without a pre-existing _k_-mer database

In the following, since no _k_-mer database was specified, qc3C will create it without prompting (`--yes`) before commencing the analysis. The resulting database will be saved to the specified output directory.
```
qc3C kmer --yes --mean-insert 300 --enzyme DpnII --reads reads_r1.fq.gz --reads reads_r2.fq.gz -output-path output
```

#### Analysing in bam mode

For BAM mode analysis, users must first create a query-name sorted BAM file of their Hi-C reads to the relevant reference sequence(s). 

Please note, that although a mapper might output reads in an apparently sorted order, this is not guaranteed. Therefore, you must <ins>explicitly sort the BAM by query-name</ins> (eg. `samtools sort -n ...`). At runtime, the BAM header is inspected for an indication of the correct order and qc3C will halt if this record is missing. Further, we encourage using BWA MEM for the mapping step as this mapper was used for development.

Steps:
1. create an index of the reference sequence
2. map the reads to reference and sort by query-name
3. run qc3C bam mode analysis
```
bwa index ref.fna.gz
bwa mem -5SP ref.fna.gz reads_r1.fq.gz reads_r2.fq.gz | samtools view -bS - | samtools sort -n -o reads2ref.bam -
qc3C bam --enzyme DpnII --fasta ref.fna.gz --bam reads2ref.bam --output-path output
```

### Results and Output files

#### Output location

All files created by qc3C are written to the location specified by the command line option `--output-path [-o]`. By default, is the current directory (`.`).

#### Console and log file

The result from each step of the analysis is directed as it completes to both the console and a log file (`qc3C.log`). The log file records all logging output and thus it is more detailed than what is sent to the console. By setting `--verbose`, the console will receive the same information as the log file.

#### Analysis report file

At the end of a run, qc3C writes a full report in both JSON (`report.qc3C.json`) and HTML (`report.qc3C.html`) formats. These files contain the same information as found in the log. While the intent of the JSON file is to ease programmatic access to the analysis results, the HTML file provides a simple structured presentation for users.

#### MultiQC support

To further aid users, we have contributed a qc3C module to [MultiQC](https://github.com/ewels/MultiQC). We encourage users to make use of this program, particularly when inspecting multiple libraries, for the benefits of its interactive and visual representation.

Both analysis modes are supported and reported separately within the MultiQC report, as well as parsing multiple experiments in to a single visual report.

**PLEASE NOTE:** MultiQC support for qc3C is currently part of the development release 1.10dev. Until such time that 1.10 is fully release, users will need to install the MultiQC master branch from github.

Using MultiQC
```
multiqc [ analysis-directory | directory-of-analysis-directories ]
```

[example reports to be added soon]

### Defining a library's digest

In either style of analysis, users must tell qc3C which enzymes were used in the DNA digestion step during library generation. This can be accomplished by either specifying these enzyme(s) by name or the commercial kit. Currently only Phase and Arima kits are defined within qc3C.

#### Using `--enzyme` (`-e`)

Enzymes are specified at runtime by their standard case-sensitive names (e.g. `DpnII`, `HindIII`, `Sau3AI`, etc). Specifying any enzyme should be possible, so long as it does not produce blunt ends. Both single and dual enzyme digests are supported, where for dual digests, users simply specify both enzymes (eg `-e Sau3AI -e MluCI`).

**Example kmer mode run for a DpnII digest and a fragment size of 300 bp**

```
qc3C kmer -m 300 -e DpnII -r reads.fq.gz
```

**Mistakes in enzyme names**

In the event a user incorrectly spells an enzyme, qc3C will attempt to suggest known ezymes with similar spelling. For instance, if a user was to accidentally use the character `1` rather than capital `I` in the name DpnII.

```
> qc3C kmer -y -m 300 -e Dpn11 -r mac_qc3C/ecoli_150_0.01.fq.gz 

INFO     | 2020-07-22 11:19:17,662 | qc3C.jellyfish | Beginning library creation
INFO     | 2020-07-22 11:19:17,663 | qc3C.jellyfish | Requested kmer size: 24
INFO     | 2020-07-22 11:19:17,663 | qc3C.jellyfish | Input FastQ files: mac_qc3C/ecoli_150_0.01.fq.gz
INFO     | 2020-07-22 11:19:17,663 | qc3C.jellyfish | Creating library: ./qc3c_kmers.jf
INFO     | 2020-07-22 11:19:22,097 | qc3C.kmer_based | Random seed was not set, using 54163895
ERROR    | 2020-07-22 11:19:22,099 |    main | Dpn11 is undefined, but its similar to: DpnII, DpnI
```

#### Using `--library-kit` (`-k`)

As a convenience, for libraries generated using either Phase or Arima kits, the enzymes can be specified as indirectly.

**Example kmer mode run generate from Phase Genomics library kit and a fragment size of 300 bp**
```
qc3C kmer -m 300 -k phase -r reads.fq.gz
```

#### Cut-site database creation

When run BAM mode, qc3C will first perform an _in silico_ digestion of the reference sequence to produce a cut-site database. For large genomes, this calculation may take a couple of minutes, but afterwards the result will be cached to a file for later reuse.

The cut-site database is saved to and retrieved from the same location as the reference sequence, therefore write access to the directory containing the reference is required.

**qc3C BAM mode run, generating the cut-site database when one does not exist**
```
qc3C bam -e DpnII -f reference.fna -b reads2ref.bam
INFO     | 2020-07-22 13:04:02,942 | qc3C.ligation | No cut-site database cache found for digest involving DpnII
INFO     | 2020-07-22 13:04:02,965 | qc3C.ligation | Building cut-site database against reference sequence and DpnII
```

### Creating a k-mer database for QC analysis

Reference-free quality assessment is a main feature of qc3C and is particularly convenient in cases where a reference is not available or of poor quality. At present, qc3C relies on [Jellyfish](https://github.com/gmarcais/jellyfish) to initially generate the k-mer database from the Hi-C read-set and afterwards uses Jellyfish's Python to analyse it.

When invoking a _k_-mer mode analysis, if you do not supply an existing _k_-mer database qc3C will offer to create it first. Otherwise, you can explicitly create the database before hand using the subcommand `mkdb` as follows.

```
qc3C mkdb -r reads_r1.fq.gz -r reads_r2.fq.gz -l kmer_db.jf 
```

Under the hood, qc3C is invoking Jellyfish to create the _k_-mer database. Jellyfish itself is quite an efficient program, but can require significant memory for larger genomes. [see Jellyfish FAQ](https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md#how-much-memory-is-needed)

### Creating A bam file for QC analysis

To quality-assess a Hi-C read-set using BAM mode, users must first map the Hi-C reads to a reference sequence obtained from the same source. 

We recommend using [BWA MEM](https://github.com/lh3/bwa) for this purpose, with options `-5SP`. Here, optns `-S` and `-P` stop BWA MEM from attempting perform mate rescue and pairing, while `-5` cause the lowest-coordinate part of split alignments to be reported as primary. This last option is particularly handy for simplifying the logic in programs analysing Hi-C read-mappings (such as bin3C). Please note that for qc3C parsing, the BAM file must be sorted by query id.
 
Below is an example of aligning Hi-C reads to a reference.  

```
bwa index ref.fa
bwa mem -5SP ref.fna.gz hic_reads.fq.gz | samtools view -bS - | samtools sort -n -o hic_to_ref.bam -
```

The state of the reference can range from a completed genome to assembly contigs, but users should keeo in mind that as the reference becomes more fragmented, assessing long-range _cis_-contacts will become increasingly hampered by short references on which to map reads. Further, _trans_-contocts (inter-molecular contacts) cannot be distinguished from those which merely connect broken segments from the same molecule.

### Suggestions on running QC analyses

#### Limit maximum observations

Analysing the entirety of a large read-set is unnecessary to obtain reliable evidence of overall quality. To limit an analysis to a smaller subset, users can specify the maximum total number of observations (`--max-obs`) and/or control acceptance rate `--sample-rate`. In testing, we have found reliable results from 100,000 reads.

#### Use an accurate insert size in _k_-mer mode

It is vitally important that an accurate estimate of the insert/fragment size is known. This value should reflect the averae size of the unknown sequence only. In testing, we have found that sequencing facilities can inadvertently quote larger values, which will affect qc3C's adjusted (or extrapolated) signal estimates.

#### Obtaining an insert size using qc3C

If a reference is available, bam mode can be used to obtain an good estimate of the median and mean observed insert size. This can subsequently be used when invoking the _k_-mer mode analysis. Using both modes may seem redudant, but each assesses your read-set from a different perspective.

While _k_-mer mode estimates the fraction of proximity-ligation pairs, while apart from estimating insert size, bam mode reports other interesting summary statistics such as the number of observed read-through events, and a breakdown of pairs (self-circle, dangling-end, etc). 

### Command line help

```$bash
usage: qc3C [-h] [-V] {bam,kmer,mkdb} ...

qc3C: Hi-C quality control

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  Version

commands:
  Valid commands

  {bam,kmer,mkdb}     choose an analysis stage for further options
```


#### Usage for _k_-mer mode

```usage: qc3C kmer [-h] [-d] [-y] [-v] [-t N] [-o PATH] [-p {range (0,1]}] [-s SEED] [-M MAX_OBS] [--no-json] [--no-html]
                 (-k {phase,arima} | -e NEB_NAME) [--ascii-base {33,64}] [--min-quality MIN_QUALITY] [--hash-size HASH_SIZE] [--kmer-size KMER_SIZE]
                 [--merged-reads] [--write-table] [--num-sample NUM_SAMPLE] [--frac-sample {range (0,1]}] [-x {range (0,1]}] -m MEAN_INSERT
                 [-l KMER_LIB] -r FASTQ_FILE

Kmer-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
  -y, --yes             Do not ask for confirmation
  -v, --verbose         Verbose output
  -t N, --threads N     Number of threads [1]
  -o PATH, --output-path PATH
                        Write output files to this folder [.]
  -p {range (0,1]}, --sample-rate {range (0,1]}
                        Sample only a proportion of all read-pairs [None]
  -s SEED, --seed SEED  Random seed used in sampling the read-set [None]
  -M MAX_OBS, --max-obs MAX_OBS
                        Stop after collecting this many observations
  --no-json             Do not write a JSON report
  --no-html             Do not write an HTML report
  -k {phase,arima}, --library-kit {phase,arima}
                        Define digest by the commercial library kit used
  -e NEB_NAME, --enzyme NEB_NAME
                        Define digest by explicitly naming up to two case-sensitive NEB enzyme names
  --ascii-base {33,64}  Ascii-encoding base for quality scores [33]
  --min-quality MIN_QUALITY
                        Minimum quality before a base position is converted to N
  --hash-size HASH_SIZE
                        Initial hash size in generating a library (eg. 10M, 2G) [10M]
  --kmer-size KMER_SIZE
                        K-mer size to use in generating a library [24]
  --merged-reads        Input reads are merged pairs
  --write-table         Save the collected observations to a file
  --num-sample NUM_SAMPLE
                        Number of samples to use in bootstrapping confidence interval [50]
  --frac-sample {range (0,1]}
                        Fraction of observations to use per-bootstrap iteration [1/3]
  -x {range (0,1]}, --max-freq-quantile {range (0,1]}
                        Ignore k-mers possessing frequencies above this quantile [0.9]
  -m MEAN_INSERT, --mean-insert MEAN_INSERT
                        Mean fragment length to use in estimating the unobserved junction rate
  -l KMER_LIB, --lib KMER_LIB
                        Jellyfish kmer database
  -r FASTQ_FILE, --reads FASTQ_FILE
                        FastQ format reads that were used to generate the k-mer library (use multiple times for multiple files)
```

#### Usage for bam mode

```usage: qc3C bam [-h] [-d] [-y] [-v] [-t N] [-o PATH] [-p {range (0,1]}] [-s SEED] [-M MAX_OBS] [--no-json] [--no-html]
                (-k {phase,arima} | -e NEB_NAME) [-q MIN_MAPQ] -f FASTA -b BAM

Alignment-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
  -y, --yes             Do not ask for confirmation
  -v, --verbose         Verbose output
  -t N, --threads N     Number of threads [1]
  -o PATH, --output-path PATH
                        Write output files to this folder [.]
  -p {range (0,1]}, --sample-rate {range (0,1]}
                        Sample only a proportion of all read-pairs [None]
  -s SEED, --seed SEED  Random seed used in sampling the read-set [None]
  -M MAX_OBS, --max-obs MAX_OBS
                        Stop after collecting this many observations
  --no-json             Do not write a JSON report
  --no-html             Do not write an HTML report
  -k {phase,arima}, --library-kit {phase,arima}
                        Define digest by the commercial library kit used
  -e NEB_NAME, --enzyme NEB_NAME
                        Define digest by explicitly naming up to two case-sensitive NEB enzyme names
  -q MIN_MAPQ, --min-mapq MIN_MAPQ
                        Minimum acceptable mapping quality [60]
  -f FASTA, --fasta FASTA
                        Reference sequences
  -b BAM, --bam BAM     Input name-sorted bam file of Hi-C reads mapped to references
```

#### Usage for mkdb
```usage: qc3C mkdb [-h] [-d] [-y] [-v] [-t N] [-o PATH] [--ascii-base {33,64}] [--min-quality MIN_QUALITY] [--hash-size HASH_SIZE]
                 [--kmer-size KMER_SIZE] -r FASTQ_FILE -l KMER_LIB

Create kmer database.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
  -y, --yes             Do not ask for confirmation
  -v, --verbose         Verbose output
  -t N, --threads N     Number of threads [1]
  -o PATH, --output-path PATH
                        Write output files to this folder [.]
  --ascii-base {33,64}  Ascii-encoding base for quality scores [33]
  --min-quality MIN_QUALITY
                        Minimum quality before a base position is converted to N
  --hash-size HASH_SIZE
                        Initial hash size in generating a library (eg. 10M, 2G) [10M]
  --kmer-size KMER_SIZE
                        K-mer size to use in generating a library [24]
  -r FASTQ_FILE, --reads FASTQ_FILE
                        FastQ format reads to use in generating the k-mer library (use multiple times for multiple files)
  -l KMER_LIB, --lib KMER_LIB
                        Output Jellyfish k-mer library base name
```
