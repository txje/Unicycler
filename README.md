# <img src="https://raw.githubusercontent.com/rrwick/Unicycler/master/misc/logo.svg" alt="logo" width="120" height="168" align="middle">Unicycler



## Table of Contents

* [Introduction](https://github.com/rrwick/Unicycler#introduction)
* [Quick usage](https://github.com/rrwick/Unicycler#quick-usage)
* [Requirements](https://github.com/rrwick/Unicycler#requirements)
* [Installation](https://github.com/rrwick/Unicycler#installation)
* [Pipeline](https://github.com/rrwick/Unicycler#pipeline)
* [Conservative, normal and bold](https://github.com/rrwick/Unicycler#conservative--normal-and-bold)
* [Options and usage](https://github.com/rrwick/Unicycler#options-and-usage)
* [Tips](https://github.com/rrwick/Unicycler#tips)
* [Unicycler align](https://github.com/rrwick/Unicycler#unicycler-align)
* [Unicycler polish](https://github.com/rrwick/Unicycler#unicycler-polish)
* [Citation](https://github.com/rrwick/Unicycler#citation)
* [Acknowledgements](https://github.com/rrwick/Unicycler#acknowledgements)
* [License](https://github.com/rrwick/Unicycler#license)



# Introduction

Unicycler is a hybrid bacterial genome assembly pipeline. By using both Illumina and PacBio/Nanopore reads, it can produce assemblies that are both accurate and complete.

Reasons to use Unicycler:
   * It has very low misassembly rates.
   * It often creates complete bacterial genomes with one circular sequence per replicon.
   * It can handle genomes with very repetitive sequences, such as _Shigella_ isolates.
   * It works with any quality of long reads - even Nanopore reads classed as 'fail' can be used as input.
   * It works with any long read depth. Approximately 10x may be required to complete a genome, but it can produce nearly-complete genomes with much fewer long reads.
   * It can be run without any long reads at all, functioning as a [SPAdes](http://bioinf.spbau.ru/spades) optimiser.
   * It produces a graph that is viewable in [Bandage](https://github.com/rrwick/Bandage), allowing for in-depth investigation.
   * It's easy to run! Unicycler is executed with only one command and no fiddling with advanced parameters is required.

Reasons to _not_ use Unicycler:
   * You only have long reads, not Illumina reads (try [Canu](https://github.com/marbl/canu) instead).
   * Your Illumina reads are poor quality (Unicycler requires a good short read assembly graph to scaffold and will not work well with low quality Illumina reads).
   * You're assembling a large eukaryotic genome or a metagenome (Unicycler is designed specifically for bacterial isolates).



# Quick usage

### Hybrid assembly
```
unicycler -1 short_reads_1.fastq.gz -1 short_reads_2.fastq.gz -l long_reads.fastq.gz -o path/to/output_dir
```

### Short read-only assembly
```
unicycler -1 short_reads_1.fastq.gz -1 short_reads_2.fastq.gz --no_long -o path/to/output_dir
```



# Requirements

* Linux or macOS
* [Python](https://www.python.org/) 3.4 or later
* A C++ compiler for installation
    * Both [GCC](https://gcc.gnu.org/) and [Clang](http://clang.llvm.org/) should work if the version isn't too old (C++11 support is required).
* [SPAdes](http://bioinf.spbau.ru/spades)

The following tools are needed for certain parts of Unicycler's pipeline. Without them, Unicycler may not be able to perform all pipeline tasks.

* [GraphMap](https://github.com/isovic/graphmap) - can accelerate long read alignment process
* [Pilon](https://github.com/broadinstitute/pilon/wiki) - required for polishing
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - required for polishing
* [Samtools](http://www.htslib.org/) version 1.0 or later - required for polishing
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) - required for rotating finished assemblies



# Installation

### Download the source code
```
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
```

### Install

If you have the necessary permissions and want to install Unicycler for all users:
```
sudo python3 setup.py install
```

If you don't have sudo access or want the installation to be specifically for your user:
```
python3 setup.py install --user
```

If that doesn't work or if you want to specify the install path:
```
python3 setup.py install --prefix=$HOME/.local
```

### Build without installation
Simply running `make` in the Unicycler directory will build the C++ code without installing Unicycler anywhere else.



# Pipeline

### Read correction

### SPAdes assembly

### Long read alignment

### Long read bridging

### Rotating completed sequences

### Short read polishing



# Conservative, normal and bold

Unicycler can be run in three modes: conservative, normal (the default) and bold. These affect the program's willingness to join sequences to make more complete assemblies.

The main difference in the modes is the quality threshold for applying bridges to the graph. In conservative mode, the threshold is set high so only the best bridges are used. This means that the assembly is less likely to be complete, but the risk of misassembly is extremely low. In bold mode, the threshold is set low so most bridges are used. This results in more complete assemblies but at an increased misassembly risk. Normal mode uses an intermediate threshold.

If the structural accuracy of the assembly is paramount to your research, conservative mode is recommended. If, however, you want a completed genome even if it contains a mistake or two, then bold mode is recommended.



# Options and usage

### Standard options

Run `unicycler --help` to view the program's most commonly used options:

```
usage: unicycler [-h] [--help_all] [--version] -1 SHORT1 -2 SHORT2 [-l LONG] [--no_long] -o OUT [--verbosity VERBOSITY] [--keep_temp KEEP_TEMP] [-t THREADS]
                 [--mode {conservative,normal,bold}] [--expected_linear_seqs EXPECTED_LINEAR_SEQS]

       __
       \ \___
        \ ___\
        //
   ____//      _    _         _                     _
 //_  //\\    | |  | |       |_|                   | |
//  \//  \\   | |  | | _ __   _   ___  _   _   ___ | |  ___  _ __
||  (O)  ||   | |  | || '_ \ | | / __|| | | | / __|| | / _ \| '__|
\\    \_ //   | |__| || | | || || (__ | |_| || (__ | ||  __/| |
 \\_____//     \____/ |_| |_||_| \___| \__, | \___||_| \___||_|
                                        __/ |
                                       |___/

Unicycler: a hybrid assembly pipeline for bacterial genomes

Help:
  -h, --help                            Show this help message and exit
  --help_all                            Show a help message with all program options
  --version                             Show Unicycler's version number

Input:
  -1 SHORT1, --short1 SHORT1            FASTQ file of short reads (first reads in each pair).
  -2 SHORT2, --short2 SHORT2            FASTQ file of short reads (second reads in each pair).
  -l LONG, --long LONG                  FASTQ or FASTA file of long reads, if all reads are available at start.
  --no_long                             Do not use any long reads (assemble with short reads only)

Output:
  -o OUT, --out OUT                     Output directory
  --verbosity VERBOSITY                 Level of stdout information (0 to 3, default: 1)
                                          0 = no stdout, 1 = basic progress indicators, 2 = extra info, 3 = debugging info
  --keep_temp KEEP_TEMP                 Level of file retention (0 to 2, default: 1)
                                          0 = only keep files at main checkpoints, 1 = keep some temp files including SAM, 2 = keep all temp files

Other:
  -t THREADS, --threads THREADS         Number of threads used to align (default: number of CPUs)
  --mode {conservative,normal,bold}     Bridging mode (default: normal)
                                          conservative = smaller contigs, lowest misassembly rate
                                          normal = moderate contig size and misassembly rate
                                          bold = longest contigs, higher misassembly rate
  --expected_linear_seqs EXPECTED_LINEAR_SEQS
                                        The expected number of linear (i.e non-circular) sequences in the underlying sequence
```

### Advanced options

Run `unicycler --helpall` to see a complete list of the program's options. These allow you to turn off parts of the pipeline, specify the location of tools (if they are not available in PATH) and adjust various settings:
```
SPAdes assembly:
  These options control the short read SPAdes assembly at the beginning of the Unicycler pipeline.

  --spades_path SPADES_PATH             Path to the SPAdes executable
  --no_spades_correct                   Skip SPAdes error correction step
  --min_kmer_frac MIN_KMER_FRAC         Lowest k-mer size for SPAdes assembly, expressed as a fraction of the read length
  --max_kmer_frac MAX_KMER_FRAC         Highest k-mer size for SPAdes assembly, expressed as a fraction of the read length
  --kmer_count KMER_COUNT               Number of k-mer steps to use in SPAdes assembly

Assembly rotation:
  These options control the rotation of completed circular sequence near the end of the Unicycler pipeline.

  --no_rotate                           Do not rotate completed replicons to start at a standard gene
  --start_genes START_GENES             FASTA file of genes for start point of rotated replicons
  --start_gene_id START_GENE_ID         The minimum required BLAST percent identity for a start gene search
  --start_gene_cov START_GENE_COV       The minimum required BLAST percent coverage for a start gene search
  --makeblastdb_path MAKEBLASTDB_PATH   Path to the makeblastdb executable
  --tblastn_path TBLASTN_PATH           Path to the tblastn executable

Pilon polishing:
  These options control the final assembly polish using Pilon at the end of the Unicycler pipeline.

  --no_pilon                            Do not use Pilon to polish the final assembly
  --bowtie2_path BOWTIE2_PATH           Path to the bowtie2 executable
  --bowtie2_build_path BOWTIE2_BUILD_PATH
                                        Path to the bowtie2_build executable
  --samtools_path SAMTOOLS_PATH         Path to the samtools executable
  --pilon_path PILON_PATH               Path to the executable Pilon Java archive file
  --java_path JAVA_PATH                 Path to the java executable
  --min_polish_size MIN_POLISH_SIZE     Sequences shorter than this value will not be polished using Pilon

Graph cleaning:
  These options control the removal of small leftover sequences after bridging is complete.

  --min_component_size MIN_COMPONENT_SIZE
                                        Unbridged graph components smaller than this size will be removed from the final graph
  --min_dead_end_size MIN_DEAD_END_SIZE
                                        Graph dead ends smaller than this size will be removed from the final graph

Long read alignment:
  These options control the alignment of long reads to the assembly graph using Graphmap and/or Unicycler-align

  --temp_dir TEMP_DIR                   Temp directory for working files ("PID" will be replaced with the process ID)
  --contamination CONTAMINATION         FASTA file of known contamination in long reads, e.g. lambda phage spike-in (default: none).
  --no_graphmap                         Do not use GraphMap as a first-pass aligner (default: GraphMap is used)
  --graphmap_path GRAPHMAP_PATH         Path to the GraphMap executable
  --scores SCORES                       Comma-delimited string of alignment scores: match, mismatch, gap open, gap extend
  --low_score LOW_SCORE                 Score threshold - alignments below this are considered poor (default: set threshold automatically)
  --min_len MIN_LEN                     Minimum alignment length (bp) - exclude alignments shorter than this length
  --keep_bad                            Include alignments in the results even if they are below the low score threshold (default: low-scoring alignments are discarded)
  --allowed_overlap ALLOWED_OVERLAP     Allow this much overlap between alignments in a single read
  --kmer KMER                           K-mer size used for seeding alignments
```



# Tips

### Running time

Unicycler is thorough and accurate, but not particularly fast. Two main factors influence the running time: the genome size/complexity and the number of long reads. Unicycler may only take an hour or so to assemble a small, simple genome with low depth long reads. On the other hand, a complex genome with many long reads may take 12 hours to finish or more. If you have a very high depth of long reads, subsampling them will make Unicycler run faster and probably not affect your assembly.

### Long read length

The length of the long reads is very important - often more so than their accuracy. A 99% accurate read isn't very useful to Unicycler if it is less than 1 kb. On the other hand, a 20 kb read is very useful, even if it is only 75% accurate. So if you do subsample your long reads, you should preferentially keep the longest ones. 

### Oxford Nanopore MinION: 1D vs 2D

Since Unicycler can tolerate low accuracy reads, 1D sequencing is probably preferable to 2D, as it will provide greater depth. However, at the time of writing, Oxford Nanopore only supports barcoding with their 2D library prep. So if you want to sequence multiple samples on a single flowcell, 2D may be the only option.

### Bad Illumina reads

Unicycler needs decent Illumina reads as input. If your Unicycler assembly looks bad, take a look at the `unbridged_graph.gfa` file in Bandage. If its in many pieces and full of dead ends, then your Illumina reads may be too poor for Unicycler to use.

### Very short contigs

Confused by very small (e.g. 1 bp) contigs? Unlike a SPAdes assembly graph where neighbouring sequences overlap by their k-mer size, a Unicycler graph has no overlaps and the sequences adjoin directly. This means that contigs in very complex parts of the graph can be quite short. They may be useless as stand-alone contigs, but they are important in the graph structure.

### Depth: chromosomes and plasmids

Unicycler normalises the depth of contigs in the graph to the median value. For a bacterial isolate, this typically means that the chromosome has a depth of around 1x and plasmids may have different (typically higher) depths.

If you see a plasmid in the assembly with a depth of 14.5x, interpret that as approximately 14 or 15 copies per cell. A plasmid with 1.3x depth means that perhaps most cells have one copy but some have more. And a plasmid with 0.7x depth might mean that it was absent from some cells in the sample.



# Unicycler align

Unicycler's algorithm for senstive semi-global alignment is available as a stand-alone alignment tool with the command `unicycler_align`.

### Semi global alignment

Semi-global alignment does not penalise end gaps so the alignment will continue until one of the two sequences ends. This can be where the two sequences overlap or where one sequence is contained within the other:

```
  TAGAA        GTGCCGGAACA         GGCCACAC     AGTAAGAT
  |||||          |||||||           |||||           |||||
ACTAGAACG        GCCGGAA       GGCTGGCCA           AAGATCTTG
```

This type of alignment is appropriate when there are no large-scale differences between the query and target sequences. For example, if you have a short read assembly graph and long reads from the same bacterial isolate (as is the case in the main Unicycler pipeline), the sequences should directly correspond and semi-global alignment is ideal.

### Compared to local alignment

Semi-global alignment may not be appropriate when the reads are being mapped to a reference genome. It can result in poor alignments around points of structural variation between the sample and the reference. For example, if the sample had a deletion relative to the reference, a read spanning that deletion would align poorly with semi-global alignment:
```
read:            AACACTAAACTTAGTCCCAA
                 |||||||||||  |   | |    
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```

In this case, local alignment (which can clip to align only part of the part) would be more appropriate:
```
read:             AACACTAAACT               TAGTCCCAA
                  |||||||||||               |||||||||
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```
[BWA-MEM](http://bio-bwa.sourceforge.net/), [LAST](http://last.cbrc.jp/) and [BLASR](https://github.com/PacificBiosciences/blasr) are all local alignment tools appropriate for such a case.



# Unicycler polish



# Citation

Paper in progress... check back later!



# Acknowledgements

Unicycler would not have been possible without [Kat Holt](https://holtlab.net/), the my fellow members of her lab and the many other researchers I work with at the University of Melbourne's [Centre for Systems Genomics](https://sysgenmelb.org/).



# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
