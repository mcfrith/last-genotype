# last-sub

This program identifies nucleotide substitutions relative to a
reference genome, from [LAST](http://last.cbrc.jp/) alignments of DNA
(or RNA) reads.

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
Now, run `last-sub` like this:

    last-sub myseq.par myseq.maf > out.txt

## Output format

The output is a table with 10 tab-delimited columns, like this:

    chr20   1118650   T   CC   9.4   1.5   CT   0    0   C~ C~ cq c~ c~
    chr20   3213195   C   AC   8.7   5.6   AT   0    0   A~ A~ C~ C~ a~ a~ c~
    chr20   3223437   G   AG   6.8   3.2   AA   13   7   A~ A~ G~ G~ a~ a~ a~ g~
    chr20   4024081   A   GG   12    2.1   AG   0    0   G! G# G3 GG g+ g. gG gM g~

* Column 1: chromosome name.

* Column 2: zero-based coordinate (i.e. the first base in the
  chromosome has coordinate 0).

* Column 3: the base at this site in the reference genome.

* Column 4: the predicted genotype.  For example, `CC` means that the
  maternal and paternal chromosomes both have `C` at this site, while
  `AC` means that one has `A` and one has `C`.

* Column 5: log10[ likelihood(predicted genotype) /
  likelihood(homozygous reference) ].

* Column 6: log10[ likelihood(predicted genotype) /
  likelihood(2nd most likely genotype) ].

* Column 7: the 2nd most likely genotype.

Columns 8 and 9 provide extra information for heterozygous sites
(i.e. where the maternal and paternal bases are predicted to differ).
For homozygous sites, they are always 0.  These columns refer to the
*phase* of this site relative to the closest preceding heterozygous
site.  In the above example, the `A` at 3213195 is predicted to lie on
the same (maternal or paternal) chromosome as the `A` at 3223437, with
the `C` at 3213195 and `G` at 3223437 on the other chromosome.

* Column 8: log10[ likelihood(predicted phase) / likelihood(the other phase) ].

* Column 9: number of reads with aligned bases at both sites.

* Column 10: the observed bases at this site, with symbols indicating
  their alignment reliability.  Each base is followed by either one
  symbol (if `-j4` was not used), or two (if `-j4` was used).  The
  first symbol, if present, comes from `-j4` and is described
  [here](http://last.cbrc.jp/doc/last-tutorial.html#example-10-ambiguity-of-alignment-columns).
  The second symbol comes from last-split and is described
  [here](http://last.cbrc.jp/doc/last-split.html#output).  Lowercase
  bases indicate reverse-strand observations; for example `g+` means
  that, actually, a `C` was observed on the reverse strand.

By default, `last-sub` only shows sites with log10[
likelihood(predicted genotype) / likelihood(homozygous reference) ] >=
6.

The output ends with a line like this:

    # Tested sites: 166730

This means that 166730 sites were covered by enough reads to have a
possibility of achieving the log likelihood threshold.

## RNA reads and pseudogenes

If your reads come from RNA, `last-sub` may make false predictions in
pseudogenes.  This is because reads from spliced genes get misaligned
to processed pseudogenes (which lack introns, allowing a misleadingly
good alignment).  You can avoid this problem by doing:

    last-sub -s50 -f1e6 myseq.par myseq.maf > out.txt

This tells it to only use reads that have typical spliced alignments:
colinear, with no splice > 10^6 bases, strong splice signals
(e.g. `GT`-`AG`), and at least one splice >= 50 bases.  (Human introns
are almost always between 50 and 10^6 bases.)

## Options

- `-h`, `--help`: show a help message and exit.

- `-m INC`, `--min=INC`: minimum increase in log10(likelihood) over
  homozygous reference (default=6).

- `-p N`, `--ploidy=N`: 1=haploid, 2=diploid, etc (default=2).

- `-f BP`, `--furthest=BP`: only use query sequences with colinear
  alignments separated by <= BP.

- `-s BP`, `--splice=BP`: only use query sequences with strong splice
  signals and least one splice >= BP.

## Gimmick

It can do arbitrary ploidy (hexaploid or whatever).

## Limitations

* Currently, `last-sub` does not use per-base quality data (e.g. from
  fastq files).  This is because it was designed for reads with more
  indel than substitution errors, such as MinION or PacBio.  I am not
  sure whether or how per-base qualities should be used for such data.

* `last-sub` does not allow for heterogeneous samples, e.g. from
  cancer, where one site may have, say, 27% `A` and 73% `C`.  But it
  will often identify such sites as heterozygous, which may be useful.

* For RNA reads: it makes no attempt to distinguish between genomic
  substitutions and RNA editing.  Also, if there is allele-biased
  expression (more RNAs from one allele), that may hide
  heterozygosity.

* It does not allow different ploidies for different chromosomes
  (e.g. X, Y, M).  This should be easy to implement, the main
  difficulty is designing the program option.

## TO DO

* Is this the best program name?
