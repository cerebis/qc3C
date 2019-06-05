# qc3C - Quality control for Hi-C DNA sequencing libraries

qc3C attempts to provide a means of assessing the proportion of "signal" within a Hi-C sequencing library. That is, read-pairs which are a true product of proximity ligation, rather than self-religation or shotgun noise. To accomplish this, two modes of analysis are available: BAM analysis and assembly-free kmer analysis. 

The BAM analysis mode requires that Hi-C reads are first mapped to a reference, preferably the same genome as was used in producing the Hi-C library. The reference can be in the form of a closed genome or assembly contigs. We recommend using [BWA MEM](https://github.com/lh3/bwa) for this purpose, with options `-5SP`.
 
As a reference genome or assembly is not necessary available at QC time, a second assembly-free analysis mode avoids this significant requirement. At present qc3C relies on a kmer database generated separately using [Jellyfish](https://github.com/gmarcais/jellyfish).

In either case, qc3C searches for evidence of proximity ligations to infer whether or not a generated library contains strong signal. Used after performing a small assessment sequencing run, this information allows researchers to choose to a full sequencing run that is more or less deep, or abandon a library entirely. 

Example usage for a library which used both Sau3AI and MluCI restriction enzymes

**BAM mode**
```$bash
> qc3C bam --mean-insert 500 -enzyme Sau3AI --enzyme MluCI --bam hic_to_ref.bam

```
**Kmer mode**
```$bash
> qc3C kmer --mean-insert 500 -enzyme Sau3AI --reads run.fq.gz --lib 24mer.jf
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