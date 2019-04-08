# qc3C - Quality control for Hi-C DNA sequencing libraries

qc3C attempts to provide a means of assessing the proportion of "signal" within a Hi-C sequencing library. That is, read-pairs which are a true product of proximity ligation, rather than self-religation or shotgun noise. To accomplish this, two modes of analysis are available: BAM analysis and assembly-free kmer analysis. 

The BAM analysis mode requires that Hi-C reads are first mapped to a reference, preferably the same genome as was used in producing the Hi-C library. The reference can be in the form of a closed genome or assembly contigs. We recommend using [BWA MEM](https://github.com/lh3/bwa) for this purpose, with options `-5SP`.
 
As a reference genome or assembly is not necessary available at QC time, a second assembly-free analysis mode avoids this significant requirement. At present qc3C relies on a kmer database generated separately using [Jellyfish](https://github.com/gmarcais/jellyfish).

In either case, qc3C searches for evidence of proximity ligations to infer whether or not a generated library contains strong signal. Used after performing a small assessment sequencing run, this information allows researchers to choose to a full sequencing run that is more or less deep, or abandon a library entirely. 

Example usage for a library which used both Sau3AI and MluCI restriction enzymes

**BAM mode**
```bash
> qc3C bam -e Sau3AI -e MluCI hic_to_ref.bam

```
**Kmer mode**
```bash
> qc3C kmer -e Sau3AI --mean-insert 500 24 run.fq.gz 24mer.jf
```


**Command help**

```$bash
usage: qc3C [-h] [-v] [-V] {bam,kmer} ...

qc3C: Hi-C quality control

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  Verbose output
  -V, --version  Version

commands:
  Valid commands

  {bam,kmer}     choose an analysis stage for further options
```

```$bash
usage: qc3C bam [-h] [--sep SEP] -e NEB_NAME [-t N] BAM

Alignment-based analysis.

positional arguments:
  BAM                   Input name-sorted bam file of Hi-C reads mapped to 
                        references

optional arguments:
  -h, --help            show this help message and exit
  --sep SEP             Delimiter to use in report table
  -e NEB_NAME, --enzyme NEB_NAME
                        Case-sensitive NEB enzyme name. Use multiple times for
                        multiple enzymes
  -t N, --threads N     Number of threads
```

```$bash
usage: qc3C kmer [-h] [--sep SEP] -e NEB_NAME [-s SEED] [-n MAX_READS] [-a]
                 [-m MEAN_INSERT] [-x MAX_COVERAGE] [-N POOL_SIZE]
                 KMER_SIZE FASTQ KMER_DB

Kmer-based analysis.

positional arguments:
  KMER_SIZE             Kmer size used in database
  FASTQ                 FastQ file used in making the kmer database
  KMER_DB               Jellyfish kmer database

optional arguments:
  -h, --help            show this help message and exit
  --sep SEP             Delimiter to use in report table
  -e NEB_NAME, --enzyme NEB_NAME
                        Case-sensitive NEB enzyme name. Use multiple times for
                        multiple enzymes
  -s SEED, --seed SEED  Random seed used in sampling the read-set
  -n MAX_READS, --max-reads MAX_READS
                        Stop after collecting N sample reads
  -a, --accept-all      Override acceptance rate and accept all useable reads
  -m MEAN_INSERT, --mean-insert MEAN_INSERT
                        Mean fragment length to use in estimating the
                        unobserved junction rate
  -x MAX_COVERAGE, --max-coverage MAX_COVERAGE
                        Ignore regions with more than this coverage
  -N POOL_SIZE, --pool-size POOL_SIZE
                        The total number of reads which are being provided for
                        consideration

```