# qc3C - reference-free quality control for Hi-C DNA sequencing libraries

The aim of qc3C is to provide a means of assessing the proportion of "signal" within a Hi-C sequencing library and thereby judge its quality as experimental data. 

By "signal" we mean the proportion of read-pairs which are true products of proximity ligation rather than pairs resulting from undesirable processes, which in downstream analyses will constitute noise in the experiment. When inferring the proportion of Hi-C pairs, no restriction is made on their separation distance. This is in contrast to reference-based methods which typically look at only the proportion of pairs that map at long-range (e.g. the number of pairs separated by more than 10 kbp). 

Although the proportion of long-range pairs is indeed an important factor, this cannot be established without the use of a reference sequence. The reference requirement is troublesome for any project which does not involve model organisms. Further, in cases where a draft assemblies is available, the degree of assembly fragmentation acts as a negative bias when trying to observe long-range pairs.

Since reference-based assessment can be of interest, qc3C can also perform this type of analysis and provide additional statistics which may be of interest. These include such read-through events, HiCPro style pair categorisation (e.g. dangling-end, self-circle, etc).

## Two modes of analysis

- **_k_-mer mode:** reference-free assessment requiring only a Hi-C read-set. 

- **BAM mode:** conventional assessment requiring read-mapping to a reference.


## Installation

### From Conda

Installation using conda is very simple.

We maintain conda packages for both qc3C and a few supporting packages, including a build of the _k_-mer counting tool Jellyfish which includes Python hooks. Please do not install the bioconda kmer-jellyfish package, as it does not possess language hooks and qc3C will fail, throwing `NoModuleFoundError: No module named 'jellyfish'`.

```$bash
conda create -n qc3c -c cerebis -c conda-forge -c bioconda qc3C
``` 

### From Docker

A docker image can be obtained and run from DockerHub.

```$bash
docker pull cerebis/qc3c:latest
docker run cerebis/qc3c:latest qc3C -h
```

The image includes binaries for the following tools:

- Jellyfish
- samtoools
- bwa
- fastp

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

#### _k_-mer mode analysis

Assu

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

### Creating A BAM file for QC analysis

To quality-assess a Hi-C read-set using BAM mode, users must first map the Hi-C reads to a reference sequence obtained from the same source. 

We recommend using [BWA MEM](https://github.com/lh3/bwa) for this purpose, with options `-5SP`. Here, optns `-S` and `-P` stop BWA MEM from attempting perform mate rescue and pairing, while `-5` cause the lowest-coordinate part of split alignments to be reported as primary. This last option is particularly handy for simplifying the logic in programs analysing Hi-C read-mappings (such as bin3C). Please note that for qc3C parsing, the BAM file must be sorted by query id.
 
Below is an example of aligning Hi-C reads to a reference.  

```$bash
bwa index ref.fa
bwa mem -5SP ref.fna.gz hic_reads.fq.gz | samtools view -bS - | samtools sort -n -o hic_to_ref.bam -
```

The state of the reference can range from a completed genome to assembly contigs, but users should keeo in mind that as the reference becomes more fragmented, assessing long-range _cis_-contacts will become increasingly hampered by short references on which to map reads. Further, _trans_-contocts (inter-molecular contacts) cannot be distinguished from those which merely connect broken segments from the same molecule.

### Creating a k-mer library for QC analysis

Reference-free quality assessment is also possible with qc3C and is convenient in cases where a reference is not available or users do not wish to spend the time creating the BAM file as above. At present, qc3C relies on [Jellyfish](https://github.com/gmarcais/jellyfish) to initially generate the k-mer library from the Hi-C read-set and afterwards uses Jellyfish's Python to analyse it.

Below is an example of creating the k-mer library from a Hi-C read-set.

```$bash
qc3C mkdb -r hic_reads.fq.gz -l hic_reads.jf 
```

### Running a QC analysis

In either mode, qc3C searches the read-pairs for evidence of proximity ligation. This evidence is then used to infer whether or not a generated library contains a strong signal (a high fraction of proximity-ligation containing read-pairs). By employing quality analysis on less expensive, small-scale pilot sequencing runs, users are better able to decide whether larger-scale sequencing should proceed and to what depth. Given sufficient budget, experiments which produced low-signal libraries might be rescued through deeper sequencing. 

Analysing the entirety of a large read-set is unnecessary to obtain reliable evidence of a read-sets overall quality. To look at only a sample of observations, users can limit the total number of observations (`--max-obs`) and/or control acceptance rate `--sample-rate`. In testing, we have found reliable results from 100,000 reads.

Example usage for a Hi-C library generated with a single enzyme.

#### BAM mode

```$bash
qc3C bam -enzyme DpnII --bam hic_to_ref.bam --fasta ref.fna.gz

```

#### _K_-mer mode

```$bash
qc3C kmer --mean-insert 500 -enzyme DpnII --reads hic_reads.fq.gz --lib hic_reads.jf
```

**Note** it is vitally important that an accurate estimate of the insert/fragment length is known. This value should purely reflect the size of the unknown sequence. Care should be taken that this value **does not** include adapter. In testing, we have found that sequencing facilities can inadvertently quote this larger value. 

If a reference is available, qc3C BAM mode will report both median and mean observed insert size, which can subsequently be used in _k_-mer mode. When available, information from both modes is of value, as _k_-mer mode estimates the fraction of proximity-ligation pairs, while apart from estimating insert size, BAM mode reports other interesting summary statistics such as the number of observed read-through events, and a breakdown of pairs (self-circle, dangling-end, etc). 

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


#### Usage for mkdb
```$bash
usage: qc3C mkdb [-h] [-d] [-v] [-t N] [-o PATH] [-s HASH_SIZE] [-k KMER_SIZE] -r FASTQ_FILE -l KMER_LIB

Create kmer database.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
  -v, --verbose         Verbose output
  -t N, --threads N     Number of threads [1]
  -o PATH, --output-path PATH
                        Write output files to this folder [.]
  -s HASH_SIZE, --hash-size HASH_SIZE
                        Initial hash size (eg. 10M, 2G) [10M]
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        Library kmer size [24]
  -r FASTQ_FILE, --reads FASTQ_FILE
                        FastQ format reads to use in generating the k-mer library (use multiple times for multiple files)
  -l KMER_LIB, --lib KMER_LIB
                        Output Jellyfish k-mer library base name
```

#### Usage for BAM mode

```$bash
usage: qc3C bam [-h] [-d] [-v] [-t N] [-o PATH] [-p {range (0,1]}] [-s SEED] [-M MAX_OBS] [--no-json] 
    [--no-html] (-k {phase,arima} | -e NEB_NAME) [-q MIN_MAPQ] -f FASTA -b BAM

Alignment-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
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

#### Usage for _K_-mer mode

```$bash
usage: qc3C bam [-h] [-d] [-v] [-t N] [-o PATH] [-p {range (0,1]}] [-s SEED] [-M MAX_OBS] [--no-json] 
    [--no-html] (-k {phase,arima} | -e NEB_NAME) [-q MIN_MAPQ] -f FASTA -b BAM

Alignment-based analysis.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Enable debug output
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
