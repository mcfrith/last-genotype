// Copyright 2017 Martin C. Frith

#include "last-genotype.hh"

#include <getopt.h>

#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

static double doubleFromString(const char *s) {
  char *e;
  double x = std::strtod(s, &e);
  if (e == s) throw std::runtime_error("bad number");
  return x;  // xxx check overflow?
}

static void run(int argc, char **argv) {
  LastGenotypeArguments args;
  args.argv = argv;
  args.min = 6;
  args.ploidy.push_back("2,chrY*:1,chrM*:1");
  args.furthest = -1;
  args.splice = -1;

  std::ostringstream help;
  help << "\
Usage: " << argv[0] << " [options] last-train.out alignments.maf\n\
\n\
Find nucleotide substitutions relative to a reference genome.\n\
\n\
Options:\n\
  -h, --help            show this help message and exit\n\
  -m INC, --min=INC     minimum increase in log10(likelihood) over homozygous\n\
                        reference (default=" << args.min << ")\n\
  -p N, --ploidy=N      1=haploid, 2=diploid, etc (default='"
       << args.ploidy[0] << "')\n\
  -f BP, --furthest=BP  only use query sequences with colinear alignments\n\
                        separated by <= BP\n\
  -s BP, --splice=BP    only use query sequences with strong splice signals\n\
                        and a splice >= BP\n\
  -V, --version         show version number and exit\n\
";

  const char sOpts[] = "hm:p:f:s:V";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "min",      required_argument, 0, 'm' },
    { "ploidy",   required_argument, 0, 'p' },
    { "furthest", required_argument, 0, 'f' },
    { "splice",   required_argument, 0, 's' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help.str();
      return;
    case 'm':
      args.min = doubleFromString(optarg);
      break;
    case 'p':
      args.ploidy.push_back(optarg);
      break;
    case 'f':
      args.furthest = doubleFromString(optarg);
      break;
    case 's':
      args.splice = doubleFromString(optarg);
      break;
    case 'V':
      std::cout << "last-genotype "
#include "version.hh"
	"\n";
      return;
    case '?':
      std::cerr << help.str();
      throw std::runtime_error("");
    }
  }

  if (optind > argc - 1) {
    std::cerr << help.str();
    throw std::runtime_error("");
  }

  args.lastTrainFile = argv[optind++];
  args.mafFiles = argv + optind;

  std::ios_base::sync_with_stdio(false);  // speed voodoo

  lastGenotype(args);
}

int main(int argc, char **argv) {
  try {
    run(argc, argv);
    if (!std::cout.flush()) throw std::runtime_error("write error");
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    const char *s = e.what();
    if (*s) std::cerr << argv[0] << ": " << s << '\n';
    return EXIT_FAILURE;
  }
}
