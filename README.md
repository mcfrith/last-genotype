# last-genotype

This program identifies nucleotide substitutions relative to a
reference genome, from [LAST](http://last.cbrc.jp/) alignments of DNA
(or RNA) reads.

# Install

You need to have standard code-building tools (e.g. git & make).  You
can get and compile `last-genotype` like this:

    git clone https://github.com/mcfrith/last-genotype.git
    cd last-genotype
    make

## Usage

First, find the substitution rates, and align the DNA (or RNA) reads,
as described
[here](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md).

It's somewhat recommended to make one modification to these alignment
recipes: add `-j4` to the `lastal` options.  `-j4` does not change the
alignments, but it adds information about the reliability of each
aligned column (and makes `lastal` slower).

In any case, you should get a file with substitution rates
(e.g. `myseq.par`), and a file of alignments (e.g. `myseq.maf`).
Now, run `last-genotype` like this:

    last-genotype -b1 myseq.par myseq.maf > out.txt

The `-b1` option tells it to require substitution evidence from both
strands.  (Omit this if e.g. your reads are from RNA forward strands
only.)

The input files may be compressed in gzip (.gz) format. You can also
pipe in the alignments, for example:

    bzcat myseq.maf.bz2 | last-genotype -b1 myseq.par > out.txt

## Output format

The output is a table with 12 tab-delimited columns, like this:

    chr20   1118650   T   CC   9.4   -0.15   CT   1.5   -0.2   0    0   C~ C~ cq c~ c~
    chr20   3213195   C   AC   8.7   -0.05   AT   5.6   0.37   0    0   A~ A~ C~ C~ a~ a~ c~
    chr20   3223437   G   AG   6.8   -0.31   AA   3.2   0.51   13   7   A~ A~ G~ G~ a~ a~ a~ g~

* Column 1: chromosome name.

* Column 2: zero-based coordinate (i.e. the first base in the
  chromosome has coordinate 0).

* Column 3: the base at this site in the reference genome.

* Column 4: the predicted genotype.  For example, `CC` means that the
  maternal and paternal chromosomes both have `C` at this site, while
  `AC` means that one has `A` and one has `C`.

* Column 5: log10[ likelihood(predicted genotype) /
  likelihood(homozygous reference) ].

* Column 6: strand bias, defined as (F - R) / (F + R), where:  
  F = L(from forward-strand alignments only)  
  R = L(from reverse-strand alignments only)  
  L = log10[ likelihood(predicted genotype) / likelihood(homozygous reference) ].  
  Note: F + R = column 5.

* Column 7: the 2nd most likely genotype.

* Column 8: log10[ likelihood(predicted genotype) /
  likelihood(2nd most likely genotype) ].

* Column 9: strand bias, defined as (F2 - R2) / (F2 + R2), where:  
  F2 = L2(from forward-strand alignments only)  
  R2 = L2(from reverse-strand alignments only)  
  L2 = log10[ likelihood(predicted genotype) / likelihood(2nd most likely genotype) ].  
  Note: F2 + R2 = column 8.

Columns 10 and 11 provide extra information for heterozygous sites
(i.e. where the maternal and paternal bases are predicted to differ).
For homozygous sites, they are always 0.  These columns refer to the
*phase* of this site relative to the closest preceding heterozygous
site.  In the above example, the `A` at 3213195 is predicted to lie on
the same (maternal or paternal) chromosome as the `A` at 3223437, with
the `C` at 3213195 and `G` at 3223437 on the other chromosome.

* Column 10: log10[ likelihood(predicted phase) / likelihood(the other
  phase) ].

* Column 11: number of reads with aligned bases at both sites.

* Column 12: the observed bases at this site, with symbols indicating
  their alignment reliability.  Each base is followed by either one
  symbol (if `-j4` was not used), or two (if `-j4` was used).  The
  first symbol, if present, comes from `-j4` and is described
  [here](http://last.cbrc.jp/doc/last-tutorial.html#example-10-ambiguity-of-alignment-columns).
  The second symbol comes from last-split and is described
  [here](http://last.cbrc.jp/doc/last-split.html#output).  Lowercase
  bases indicate reverse-strand observations; for example `g+` means
  that, actually, a `C` was observed on the reverse strand.

By default, `last-genotype` only shows sites with
log10[ likelihood(predicted genotype) / likelihood(homozygous reference) ] >=
6.

The output ends with a line like this:

    # Tested sites: 166730

This means that 166730 sites were covered by enough reads to have a
possibility of achieving the log likelihood threshold.

## RNA reads and pseudogenes

If your reads come from RNA, `last-genotype` may make false
predictions in pseudogenes.  This is because reads from spliced genes
get misaligned to processed pseudogenes (which lack introns, allowing
a misleadingly good alignment).  You can avoid this problem by doing:

    last-genotype -s50 -f1e6 myseq.par myseq.maf > out.txt

This tells it to only use reads that have typical spliced alignments:
colinear, with no splice > 10^6 bases, strong splice signals
(e.g. `GT`-`AG`), and at least one splice >= 50 bases.  (Human introns
are almost always between 50 and 10^6 bases.)

## Big input data

`last-genotype` needs to hold intermediate data of comparable size to
the input.  If this is very large, temporary files will be used.  It's
necessary to put the temporary files on a large enough (and preferably
fast) disk: you can control this with the `-T` option or `TMPDIR`
environment variable.  The temporary files are not visible as normal
files, but you can measure disk usage with a command such as `df -h`.

## Options

- `-h`, `--help`: show a help message and exit.

- `-m INC`, `--min=INC`: minimum increase in log10(likelihood) over
  homozygous reference (default=6).

- `-b INC`, `--bias=INC`: require that the column-6 strand bias has
  magnitude < BIAS.

- `-p N`, `--ploidy=N`: 1=haploid, 2=diploid, etc
  (default=`'2,chrY*:1,chrM*:1'`).

- `-f BP`, `--furthest=BP`: only use query sequences with colinear
  alignments separated by <= BP.

- `-s BP`, `--splice=BP`: only use query sequences with strong splice
  signals and least one splice >= BP.

- `-c FILE`, `--class-file=FILE`: write a classification of the query
  sequences by maternal/paternal chromosome.

- `-S SIZE`, `--buffer-size=SIZE`: maximum amount of memory to use
  (roughly).  E.g. `32G` means 32 GibiBytes.

- `-T DIR`, `--temporary-directory=DIR`: put temporary files in this
  directory.  The default is to use the `TMPDIR` environment variable,
  or if that is not specified, `/tmp`.

- `-v`, `--verbose`: show progress messages.

## Ploidy

You can specify ploidies with option `-p`.  Examples:

* `-p1` : all chromosomes are haploid, i.e. just 1 copy.

* `-p2` : all chromosomes are diploid, i.e. 2 copies (maternal and paternal).

* `-p6` : all chromosomes are hexaploid.

* `-p'chr21:3'` : `chr21` is triploid.

* `-p'chr21*:3'` : all chromosomes whose names start with `chr21`
  (e.g. `chr21_GL383579v2_alt`) are triploid.

* `-p'2,chrW*:1,chrM*:1'` : all chromosomes are diploid, except those
  whose names start with `chrW` or `chrM`, which are haploid.

This option accepts one or more comma-separated ploidy specifications.
Later specifications override earlier ones.  Your specifications are
*appended* to the default ones.

You can specify chromosome names using these symbols:

    *        zero or more of any character
    ?        any single character
    [abc]    any character in abc
    [!abc]   any character not in abc

## Classifying the reads

If you use the class-file option:

    last-genotype -c classify.txt -b1 myseq.par myseq.maf > out.txt

`last-genotype` will write a file showing which (maternal or paternal)
chromosome each read is from.  This file has 10 tab-delimited columns:

    readJ   799     1453    chr11   838132  838817  +       a       12      7
    readK   40      119     chr11   838728  838817  +       b       0.89    1
    readL   70      370     chr11   838491  838819  -       b       7.1     4

Each line corresponds to one of the input alignments.  (So it is
possible for one read to appear in more than one line, if it was split
into more than one alignment.)  The first 3 columns show the aligned
part of the read, in [BED3
format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).  The next
3 columns show the aligned part of the reference genome.  The 7th
column indicates the DNA strand.

The 8th column shows which (maternal or paternal) chromosome the read
is from.  `a` means "left" and `b` means "right".  For example, given:

    chr11   838721  T       TC      ...
    chr11   838759  T       AT      ...

`a` indicates the chromosome that has `T` at 838721 and `A` at 838759,
whereas `b` indicates the chromosome that has `C` at 838721 and `T` at
838759.

The 9th column is: log10[ likelihood(predicted chromosome) /
likelihood(the other chromosome) ].

The 10th column shows the number of bases aligned to heterozygous
sites.

Cases with zero in column 9 or 10 are not written, because they cannot
be classified.

Blank lines are written at locations with unknown phasing (zero in
column 10 of the main `last-genotype` output).  Thus, at each blank
line, there is a 50/50 chance that `a` and `b` (i.e. left and right)
switch to the opposite (maternal or paternal) chromosome.  In other
words, the classification is meaningful only within each
blank-line-separated block of lines.

## Limitations

* It only finds substitutions (not deletions, insertions, inversions,
  etc.)

* Currently, `last-genotype` does not use per-base quality data
  (e.g. from fastq files).  This is because it was designed for reads
  with more indel than substitution errors, such as MinION or PacBio.
  I am not sure whether or how per-base qualities should be used for
  such data.

* `last-genotype` does not allow for heterogeneous samples, e.g. from
  cancer, where one site may have, say, 27% `A` and 73% `C`.  But it
  will often identify such sites as heterozygous, which may be useful.

* Phasing (and classifying reads) is attempted only for diploid
  chromosomes.

* For RNA reads: it makes no attempt to distinguish between genomic
  substitutions and RNA editing.  Also, if there is allele-biased
  expression (more RNAs from one allele), that may hide
  heterozygosity.

* It doesn't allow varying ploidy *within* a chromosome
  (e.g. pseudo-autosomal regions).

* It can't combine datasets with different substitution rates.
