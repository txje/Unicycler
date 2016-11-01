<p align="center"><img src="misc/logo.png" alt="Unicycler" width="600" height="210"></p>

Unicycler is a hybrid assembly pipeline for bacterial genomes. It uses both Illumina and PacBio/Nanopore reads to produce complete and accurate assemblies.



# Table of contents

* [Introduction](#introduction)
* [Quick usage](#quick-usage)
* [Requirements](#requirements)
* [Installation](#installation)
* [Pipeline](#pipeline)
* [Conservative, normal and bold](#conservative-normal-and-bold)
* [Options and usage](#options-and-usage)
* [Output files](#output-files)
* [Tips](#tips)
* [Unicycler align](#unicycler-align)
* [Unicycler polish](#unicycler-polish)
* [Citation](#citation)
* [Acknowledgements](#acknowledgements)
* [License](#license)



# Introduction

Reasons to use Unicycler:
   * It has very low misassembly rates.
   * It often completes bacterial genomes with one circular sequence per replicon.
   * It handles organisms with highly repetitive genomes, such as _Shigella_.
   * It works with long reads of any quality - even Nanopore reads classed as 'fail' can be used as input.
   * It works with any long read depth. Approximately 10x may be required to complete a genome, but it can make nearly-complete genomes with far fewer long reads.
   * It can be run without any long reads, functioning as a [SPAdes](http://bioinf.spbau.ru/spades) optimiser.
   * It produces assembly graphs, allowing for in-depth investigations in [Bandage](https://github.com/rrwick/Bandage).
   * It's easy to use, runs with just one command and doesn't require tinkering with parameters.

Reasons to __not__ use Unicycler:
   * You only have long reads, not Illumina reads (try [Canu](https://github.com/marbl/canu) instead).
   * Your Illumina reads are poor quality (Unicycler requires a good short read assembly graph - [more info here](#bad-illumina-reads)).
   * You're assembling a large eukaryotic genome or a metagenome (Unicycler is designed for bacterial isolates).
   * You're very impatient (Unicycler is not as fast as alternative assemblers).



# Quick usage

### Hybrid assembly
```
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads.fastq.gz -o path/to/output_dir
```

### Short read-only assembly
```
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --no_long -o path/to/output_dir
```



# Requirements

* Linux or macOS
* [Python](https://www.python.org/) 3.4 or later
* C++ compiler
    * [GCC](https://gcc.gnu.org/), [Clang](http://clang.llvm.org/) and [ICC](https://software.intel.com/en-us/c-compilers) should all work if the version isn't too old (C++11 support is required).
* [SPAdes](http://bioinf.spbau.ru/spades)

Unicycler needs the following tools for certain parts of its pipeline. They are optional, but without them Unicycler will not be able to perform all pipeline tasks:

* [GraphMap](https://github.com/isovic/graphmap) - can accelerate long read alignment process
* [Pilon](https://github.com/broadinstitute/pilon/wiki) - required for polishing
* Java - required for polishing
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - required for polishing
* [Samtools](http://www.htslib.org/) version 1.0 or later - required for polishing
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) - required for rotating finished assemblies



# Installation

### Typical installation
```
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install
```
If the last command complains about permissions, you may need to run it with `sudo`.

### Other installation commands

* Install just for your user: `python3 setup.py install --user`
    * If you get a strange 'can't combine user with prefix' error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Unicycler`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Unicycler.git`
* Install with specific Makefile options: `python3 setup.py install --makeargs "CXX=icpc"`
* Build and run without installing:
    * `make`
    * `./unicycler-runner.py`

# Pipeline

### Read correction

Unicycler uses SPAdes to perform read error correction before assembling the Illumina reads. This can be disabled with `--no_correct` if your Illumina reads are very high quality or you've already performed read QC.

### SPAdes assembly

Unicycler uses SPAdes to assembly the Illumina reads into an assembly graph. It tries assemblies at a wide range of k-mer sizes, evaluating the graph at each one. It chooses the assembly graph which best balances contig length and dead end count. If the Illumina reads are high quality, this will result in an assembly graph with long contigs but few to no dead ends.

Unicycler also cleans the SPAdes assembly graphs, filtering out contigs which are very low depth and likely to be due to contamination or errors.

### Determine multiplicity

In future steps, Unicycler will scaffold the graph using SPAdes contigs and long reads. To do this, it must distinguish between single-copy contigs and collapsed repeats. It does this with a greedy algorithm that takes both read depth and graph connectivity. This process finds single-copy contigs not only in the bacterial chromosome but also in plasmids of any read depth.

### SPAdes bridging

At this point, the assembly graph does not contain the SPAdes repeat resolution. To apply this to the graph, Unicycler builds bridges between single-copy contigs using the information in the SPAdes `contigs.paths` file. These are applied to the graph to make the `spades_bridges_applied.gfa` output - the most resolved graph Unicycler can make using only the Illumina reads.

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

Run `unicycler --help_all` to see a complete list of the program's options. These allow you to turn off parts of the pipeline, specify the location of tools (only necessary if they are not in PATH) and adjust various settings:
```
SPAdes assembly:
  These options control the short read SPAdes assembly at the beginning of the Unicycler pipeline.

  --spades_path SPADES_PATH             Path to the SPAdes executable
  --no_correct                          Skip SPAdes error correction step
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



# Output files

Depending on the input files and the value used for `--keep_temp`, Unicycler may only only produce some of these. Also, all outputs except for assembly.gfa and assembly.fasta will be prefixed with a number so they are in chronological order.

File                           | Description
------------------------------ | ---------------------------------------------------------------------------
unbridged_graph.gfa            | short read assembly graph before any bridges have been applied
spades_bridges_applied.gfa     | SPAdes bridges applied, before any cleaning or merging
cleaned.gfa                    | redundant contigs removed from the graph
merged.gfa                     | contigs merged together where possible
long_read_bridges_applied.gfa  | Long read bridges applied, before any cleaning or merging
cleaned.gfa                    | redundant contigs removed from the graph
merged.gfa                     | contigs merged together where possible
final_clean.gfa                | more redundant contigs removed
rotated.gfa                    | circular replicons rotated and/or flipped to a start position
polished.gfa                   | after a round of Pilon polishing
assembly.gfa                   | __the final assembly in graph format__
assembly.fasta                 | __the final assembly in FASTA format__ (exact same contigs as assembly.gfa)



# Tips

### Running time

Unicycler is thorough and accurate, but not particularly fast. Two main factors influence the running time: the genome size/complexity and the number of long reads. Unicycler may only take an hour or so to assemble a small, simple genome with low depth long reads. On the other hand, a complex genome with many long reads may take 12 hours to finish or more. If you have a very high depth of long reads, subsampling them will make Unicycler run faster and probably not affect your assembly.


### Long read length

The length of the long reads is very important - often more so than their accuracy. A 99% accurate read isn't very useful to Unicycler if it is less than 1 kb. On the other hand, a 20 kb read is very useful, even if it is only 75% accurate. So if you do subsample your long reads, you should preferentially keep the longest ones.


### Poretools

[Poretools](http://poretools.readthedocs.io/en/latest/) can turn your Oxford Nanopore FAST5 reads into a FASTQ file appropriate for Unicycler. Here's an example command:
```
poretools fastq --type best --min-length 1000 path/to/fast5/dir/ > nanopore_reads.fastq
```
If you have 2D reads, the `--type best` option makes Poretools give only one FASTQ read per FAST5 file (if you have 1D reads, you can exclude that option). Adjust the `--min-length 1000` parameter to suit your dataset - a larger value would be appropriate if you have lots of long reads.


### Oxford Nanopore: 1D vs 2D

Since Unicycler can tolerate low accuracy reads, [Oxford Nanopore 1D sequencing](https://nanoporetech.com/applications/dna-nanopore-sequencing) is probably preferable to 2D, as it can provide twice as many reads. However, at the time of writing, the 2D library prep supports barcoding. So if you want to sequence multiple samples on a single flowcell, 2D is currently the only option.


### Bad Illumina reads

Unicycler needs decent Illumina reads as input - ideally with uniform read depth and no unsequenced parts.

You can look at the `unbridged_graph.gfa` file (the first graph Unicycler saves to file) in Bandage to get a quick impression of the Illumina read quality:

<p align="center"><img src="misc/illumina_graph_comparison.png" alt="Graphs of varying quality"></p>

__A__ is an very good Illumina read graph - the contigs are long and there are no dead ends. This read set is ideally suited for use in Unicycler.

__B__ is also a good graph. The genome is more complex, resulting in a more tangled structure, but there are still very few dead ends (you can see one in the lower left). This read set would also work well in Unicycler.

__C__ is a disaster! It is broken into many pieces, probably because parts of the genome got no read depth at all. While you can still use Unicycler to resolve this assembly with long reads, the risk of small errors and misassemblies is considerably higher.


### Very short contigs

Are you confused by very small (e.g. 2 bp) contigs in Unicycler assemblies? Unlike a SPAdes graph where neighbouring sequences overlap by their k-mer size, a Unicycler graph has no overlaps and the sequences adjoin directly. This means that contigs in very complex parts of the graph can be quite short. They may be useless as stand-alone contigs but are still very important in the graph structure.

<p align="center"><img src="misc/short_contigs.png" alt="Short contigs in assembly graph"></p>

If you're curious, this example is the rDNA region of a bacterial genome assembly graph. There are seven rDNA copies with regions that are the same (assembly collapsed these into single contigs) but in some places differ (leading to divergences in the graph like this one).


### Depth: chromosomes and plasmids

Unicycler normalises the depth of contigs in the graph to the median value. For a bacterial isolate, this typically means that the chromosome has a depth of around 1x and plasmids may have different (typically higher) depths.

(IMAGE OF VARIOUS PLASMID DEPTHS)

In this example, plasmid __A__ occurs in approximately 14 or 15 copies per cell in the sample. Most cells probably had one copy of plasmid __B__ but some had more. Since sequencing biases (such as GC bias) can affect read depth, these relative depths should be interpreted loosely.



# Unicycler align

Unicycler's algorithm for senstive semi-global alignment is available as a stand-alone alignment tool with the command `unicycler_align`.
```
unicycler_align --reads queries.fastq --ref target.fasta --sam output.sam
```

### Semi-global alignment

Semi-global alignment (a.k.a. glocal, overlap or free end-gap alignment) will not clip an alignment until one of the two sequences ends. This can be where one sequence is contained within the other or where the two sequences overlap:
```
  TAGAA        GTGCCGGAACA         GGCCACAC     AGTAAGAT
  |||||          |||||||           |||||           |||||
ACTAGAACG        GCCGGAA       GGCTGGCCA           AAGATCTTG
```

In contrast, local alignment will align only the best matching parts, clipping the alignment where the quality becomes poor:
```
      CGAACAGCATACTTG
          ||||||||
ACGTCAGACTCAGCATACGCATCTAGA
```

Semi-global alignment is appropriate when there are no structural differences between the query and reference sequences. For example, when you have a short read assembly graph and long reads from the same bacterial isolate (as is the case in the Unicycler pipeline). In this scenario, there may be small scale differences (due to read errors) but no large scale differences, and semi-global alignment is ideal.

### Compared to local alignment

Semi-global alignment is probably not appropriate for mapping reads to a more distant reference genome. It does not cope with points of structural variation between the sample and the reference. For example, if the sample had a deletion relative to the reference, a read spanning that deletion would align poorly with semi-global alignment:
```
read:            AACACTAAACTTAGTCCCAA
                 |||||||||||  |   | |    
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```

Local alignment (which can align only part of the read) would be more appropriate:
```
read:            AACACTAAACT               TAGTCCCAA
                 |||||||||||               |||||||||
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```
Try [BWA-MEM](http://bio-bwa.sourceforge.net/), [LAST](http://last.cbrc.jp/) or [BLASR](https://github.com/PacificBiosciences/blasr) if you need a local alignment tool.



# Unicycler polish

Unicycler polish is a script to repeatedly polish a completed assembly using all available reads. It can be given Illumina reads, long reads or (ideally) both. When both Illumina and long reads are available, Unicycler polish can fix assembly errors, even in repetitive parts of the genome which cannot be polished by short reads alone.

Unicycler polish uses an exhaustive iterative process that is time-consuming but can be necessary to resolve the sequence in repeat regions. For example, consider a genome with two very similar regions, A and B, and there are assembly errors in both. Polishing is initially difficult because the errors may cause reads which should map to A to instead map to B and vice versa. However, after some of these errors are fixed, more reads will map to their correct locations, allowing for more errors to be fixes, allowing more reads to map correctly, etc.

### Process
1. If Illumina reads are available:
    1. Run [Pilon](https://github.com/broadinstitute/pilon/wiki) in 'bases' mode (substitutions and small indels). If any changes were suggested, apply them and repeat this step.
    2. Run Pilon in 'local' mode (larger variants), and assess each change with ALE. If any variant improves the ALE score, apply it and go back to step 1-i.
2. If long reads are available:
    1. Run [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)/[Nanopolish](https://github.com/jts/nanopolish) and gather all suggested small changes.
    2. Use [FreeBayes](https://github.com/ekg/freebayes) to assess each long read-suggested change by looking for ambiguity in the Illumina read mapping. If any were found, apply them and go back to step 2-i.
3. If Illumina reads are available:
    1. Execute step 1 again.
    2. Run Pilon/GenomicConsensus/Nanopolish again (all that apply) and assess each suggested variant with ALE. If any improves the ALE score, apply it and repeat this step.

### Example commands
Polishing with only Illumina reads:
```
unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -a assembly.fasta
```
Polishing with only PacBio reads:
```
unicycler_polish --pb_bam subreads.bam -a assembly.fasta
unicycler_polish --pb_bax path/to/*bax.h5 -a assembly.fasta
```
Polishing with both Illumina and PacBio reads:
```
unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --pb_bam subreads.bam -a assembly.fasta
unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --pb_bax path/to/*bax.h5 -a assembly.fasta
```



# Citation

Paper in progress... check back later!



# Acknowledgements

Unicycler would not have been possible without [Kat Holt](https://holtlab.net/), my fellow researchers in her lab and the many other people I work with at the University of Melbourne's [Centre for Systems Genomics](https://sysgenmelb.org/).



# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
