# qc3C - reference-free quality control for Hi-C DNA sequencing libraries

The primary aim of qc3C is to infer the proportion of "signal" within a Hi-C sequencing library and thereby provide a means of quality control. A big advantange of qc3C is that this can be done without access to a reference sequence, which until now has been a significant stopping point for projects not involving model organisms. 

By "signal" we mean the proportion of read-pairs which are true products of proximity ligation rather than pairs resulting from undesirable processes, which in downstream analyses will constitute noise in the experiment. When inferring the proportion of Hi-C pairs, no restriction is made on their separation distance. This is in contrast to reference-based methods which typically look at only the proportion of pairs that map at long-range (e.g. the number of pairs separated by more than 10 kbp). 

Although the proportion of long-range pairs is indeed an important attribute, this type of information cannot be established without the use of a reference sequence (preferably high-quality). In cases where a draft assembly is used, sequence fragmentation reduces the likelihood of observing pairs with large separation.

Still, details from reference-based assessment remain of interest and therefore qc3C can also perform this type of analysis. What we call "bam mode" analysis. Statistics obtained from this mode include such details as the number of read-through events and HiCPro style pair categorisation (e.g. dangling-end, self-circle, etc).

## Two modes of analysis

- **_k_-mer mode:** reference-free assessment requiring only a Hi-C read-set. 

- **bam mode:** conventional assessment requiring read-mapping to a reference.

## Installation

### Requirements

Due to our dependency on Jellyfish, qc3C requires Linux x86_64 to run natively. With that said, we maintain a Docker image which will run on other platforms.

### Using Conda

We maintain conda packages for both qc3C and a few supporting packages on the Anaconda channel [cerebis](https://anaconda.org/cerebis). 

Installation is accomplished as follows.
```
conda create -y -n qc3c -c cerebis -c conda-forge -c bioconda qc3C
``` 

Note: our anaconda channel includes a build of the _k_-mer counting tool Jellyfish with Python language hooks. **Please do not install the bioconda kmer-jellyfish package** as it does not possess language hooks and qc3C will fail to run, throwing `NoModuleFoundError: No module named 'jellyfish'`.

### Using Docker

A docker image is available from DockerHub at docker://cerebis/qc3c.

```
docker pull cerebis/qc3c
docker run cerebis/qc3c -h
```

### Using Singularity

We rely on singularity users importing our Docker image, which could be done as follows.

```
singularity build qc3C.sif docker://cerebis/qc3c
singularity run qc3C.sif
```

### From Github

qc3C can be installed directly from github, however this requires that [Jellyfish](https://github.com/gmarcais/Jellyfish), along with its Python hooks be installed first. Although the basic Jellyfish binaries are easily built and installed, a build which correctly creates the Python hooks is more problematic. To remedy this, users are encouraged to use our Jellyfish conda package, after which qc3C is easily installed using Pip.

**Note:** Do not use Bioconda's Jellyfish package, as it contains only the Jellyfish binaries and no language hooks. As a result, qc3C will fail to run with the error `NoModuleFoundError: No module named 'jellyfish'`. 

The following bash snippet prepares the conda environment with our Jellyfish package and installs qc3C directly from Github.

```$bash
conda create -y -n qc3c -c cerebis -c conda-forge -c bioconda jellyfish
conda activate qc3c
pip install git+https://github.com/cerebis/qc3C
```

## Using qc3C

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

Users must first map reads to the chosen reference sequence and create a query-name ordered bam file. We strongly encourage the use bwa mem for the step of mapping. _Note, qc3C currently will not create the bam file for you._

1. first prior to mapping, an index of the reference sequence is created
2. next, reads are mapped to reference using bwa mem and sorted by query-name
3. last, bam mode analysis is run
```
bwa index ref.fna.gz
bwa mem -5SP ref.fna.gz reads_r1.fq.gz reads_r2.fq.gz | samtools view -bS - | samtools sort -n -o reads2ref.bam -
qc3C bam --enzyme DpnII --fasta ref.fna.gz --bam reads2ref.bam --output-path output
```

#### Results

In each case within the output directory (default is `.`), qc3C will write the analysis results to the console as well as to the log file `qc3C.log`. In addition, a report is written in both JSON and HTML formats, where the HTML report is a very simple conversion of the JSON report to a user-friendly layout, while the JSON report is intended for machine parsing.

To further aid users, we have contributed to the visual report tool [MultiQC](https://github.com/ewels/MultiQC), so that qc3C reports can be inspected in a more visual format. The qc3C reports are automatically recognised by MultiQC.

Using the standard interface of MultiQC, in the output folder as above
```
multiqc output
```

### Defining a library's digest

In either style of analysis, users must tell qc3C which enzymes were used in the DNA digestion step during library generation. This can be accomplished by either specifying these enzyme(s) by name or the commercial kit. Currently only Phase and Arima kits are defined within qc3C.

#### --enzyme (-e)

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

#### --library-kit (-k)

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


### Notes on running QC analyses

#### Limit maximum observations

Analysing the entirety of a large read-set is unnecessary to obtain reliable evidence of overall quality. To limit an analysis to a smaller subset, users can specify the maximum total number of observations (`--max-obs`) and/or control acceptance rate `--sample-rate`. In testing, we have found reliable results from 100,000 reads.

#### Accurate insert size in _k_-mer mode

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
