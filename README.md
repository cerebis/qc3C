# qc3C - Quality control for Hi-C libraries

qc3C analyzes Hi-C reads for evidence of proximity ligations to infer whether tha dataset contains strong signal.

As input, qc3C takes a BAM file of mapped Hi-C read-pairs to the desired reference. The reference may be any type of sequence, from closed genomes to metagenomic assembly contigs.

The analysis currently involves inspecting the mapped reads for evidence of cut-site terminated alignments, proximity ligation junctions and split alignments. These features are strong indicators of the occurence of Hi-C proximity ligation events, but do not inform us about the distribution of magnitude of separation between pairs. Therefore a second stage of analysis is still warranted which assesses the range of separation, potentially ratio of trans to cis ligations, and spurious ligation products (noise).

Most enzymes should be known to the software by name (case-sensitive).

Example usage for a library which used both Sau3AI and MluCI restriction enzymes
```bash
> qc3C -e Sau3AI -e MluCI hic_to_ref.bam
```

Help
```
usage: qc3C [-h] [--log] [-t THREADS] -e NEB_NAME BAM

positional arguments:
  BAM                   Input bam file of Hi-C reads mapped to references

optional arguments:
  -h, --help            show this help message and exit
  --log                 Keep a log of pairs with cut-sites
  -t N, --threads N     Number of threads
  -e NEB_NAME, --enzyme NEB_NAME
                        Case-sensitive NEB enzyme name. Use multiple times for
                        multiple enzymes
```
