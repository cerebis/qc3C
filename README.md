# qc3C - quality control for Hi-C DNA sequencing libraries

The aim of qc3C is to provide a means of assessing the proportion of "signal" within a Hi-C sequencing library and thereby judge its quality as experimental data. 

We use "signal" to refer to the proportion of read-pairs which are true products of proximity ligation, rather than pairs that result from undesirable processes that in downstream analysis will constitute the noise in the experiment. 
 
To accomplish this, two modes of analysis are available:

- BAM mode

    Conventional assessment requiring a reference sequence and read-mapping.

- _K_-mer mode

    **Reference-free** assessment requiring only a Hi-C read-set. 

## Installation

### From Conda

Installation using conda is very simple.

We maintain conda packages for both qc3C and a few supporting packages, including a build of the _k_-mer counting tool Jellyfish which includes Python hooks. Please do not install the bioconda kmer-jellyfish package, as this does not possess any language hooks and qc3C will throw `NoModuleFoundError: No module named 'jellyfish'`.

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

qc3C can be installed directly from github, however this requires that [Jellyfish](https://github.com/gmarcais/Jellyfish), along with its Python hooks be installed first. Although the basic Jellyfish binaries are easily built and installed, the a build which correctly creates the Python hooks is more problematic. To remedy this, users are encouraged to use our Jellyfish conda package, after which qc3C is easily installed using Pip.

**Note:** Do not use Bioconda's Jellyfish package, as it contains only the Jellyfish binaries and no language hooks. As a result, when run qc3C will throw `NoModuleFoundError: No module named 'jellyfish'`. 

The following bash snippet prepares the conda environment with our Jellyfish package and installs qc3C directly from Github.

```$bash
conda create -y -n qc3c -c cerebis -c bioconda -c conda-forge jellyfish
conda activate qc3c
pip install git+https://github.com/cerebis/qc3C
```

## Using qc3C

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
