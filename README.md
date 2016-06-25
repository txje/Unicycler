# Unicycler

This repository contains three separate (but related) tools:
* Unicycler: a bacterial genome assembler for hybrid read sets
* Antpath: a sensitive, semi-global long read aligner
* Scrutinate: an assembly checker



# Installation

Unicycler, Antpath and Scrutinate are all built/installed together. Simply run: `python3 setup.py install`.



# Unicycler

Unicycler is a bacterial genome assembler for use with hybrid read sets (both short and long reads). It uses [SPAdes](http://bioinf.spbau.ru/spades) to construct an assembly graph from short reads, then it uses long read alignments to scaffold and simplify the graph.

### Input requirements

Unicycler requires short reads (e.g. Illumina read) of decent sequencing depth. Without a good set of short reads, SPAdes will not be able to make a reasonably complete assembly graph.

The long reads, however, can be of any depth. Moderate long read depth should be sufficient to complete a bacterial genome assembly. With very low long read depth, Unicycler may not be able to complete an assembly, but it can still use them to improve the assembly.



# Antpath

Antpath sensitively aligns error-prone long reads (e.g. PacBio or Nanopore) to one or more references in a semi-global manner.

Semi-global alignment allows for unpenalised end gaps, but the alignment will continue until one of the two sequences ends. This includes cases where the two sequences overlap and cases where one sequence is contained within the other:

```
  AAAAA        AAAAAAAAAAA         AAAAAAAA     AAAAAAAA
  |||||          |||||||           |||||           |||||
BBBBBBBBB        BBBBBBB       BBBBBBBBB           BBBBBBBBB
```

Antpath is intended for cases where the reads and reference are expected to match perfectly (or at least as perfectly as error-prone long reads can match). An example of an appropriate case would be if the reference sequences are assembled contigs of a bacterial strain and the long reads are from the same strain.

Required inputs:
  1) FASTA file of one or more reference sequences
  2) FASTQ file of long reads

Output: SAM file of alignments

### Required arguments:
* `--ref`: FASTA file containing one or more reference sequences
* `--reads`: FASTQ file of long reads
* `--sam`: SAM file of resulting alignments

### Some important optional arguments:
* `--no_graphmap`: If you use this flag, the aligner will align all reads by itself. If you don't use this flag (i.e. the default behaviour), the aligner will first align reads using GraphMap, which is faster. Then it will only try to align reads for which GraphMap failed or overlap the end of a reference.
* `--graphmap_path`: Use this argument to specify the location of the GraphMap executable. If not used, this script will just try `graphmap` (i.e. assumes that GraphMap is in your PATH).
* `--scores`: Comma-delimited string of alignment scores: 'match, mismatch, gap open, gap extend'. Examples: '3,-6,-5,-2' (the default), '5,-6,-10,0' (BLASR's default), '2,-5,-2,-1' (BWA-MEM's default), '2,-3,-5,-2' (blastn's default).
* `--min_len`: Alignments shorter than this length (in bp) will not be included in the results.
* `--threads`: The number of CPU threads to use (default is to use all available CPUs).
* `--verbosity`: How much stdout to produce. 0 = no stdout. 1 = default. 2 = extra output (display the alignments for each read). 3 = lots of extra output (display the attempted alignments for each read). 4 = tons (I use this for debugging - you probably shouldn't use it).



# Scrutinate

Scrutinate allows a user to assess whether a completed bacterial genome assembly is correct, using long reads. It requires long read alignments (produced by `semi_global_long_read_aligner.py`) and produces summary tables and plots.
