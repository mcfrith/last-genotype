// Copyright 2017 Martin C. Frith

#include "last-genotype.hh"

#include <getopt.h>

#include <cfloat>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

const size_t maxSize = -1;

static double doubleFromString(const char *s) {
  char *e;
  double x = std::strtod(s, &e);
  if (e == s) throw std::runtime_error("bad number");
  return x;  // xxx check overflow?
}

static size_t sizeFromString(const std::string &s) {
  std::istringstream iss(s);
  size_t x;
  iss >> x;
  if (!iss) throw std::runtime_error("bad size: " + s);
  char suffix = 'K';
  iss.get(suffix);
  const char *i = "bKMGTP";
  const char *j = std::strchr(i, suffix);
  if (!j || iss.get(suffix)) throw std::runtime_error("bad size: " + s);
  for (; i < j; ++i) {
    if (x > maxSize / 1024) throw std::runtime_error("too big: " + s);
    x *= 1024;
  }
  return x;
}

static size_t defaultBufferSize() {
  size_t s = 16;
  for (int i = 0; i < 30 && s <= maxSize / 2; ++i) s *= 2;
  return s;
}

static void run(int argc, char **argv) {
  LastGenotypeArguments args;
  args.argv = argv;
  args.min = 6;
  args.bias = DBL_MAX;
  args.ploidy.push_back("2,chrY*:1,chrM*:1");
  args.furthest = -1;
  args.splice = -1;
  args.buffer_size = defaultBufferSize();
  args.temporary_directory = 0;
  args.verbose = 0;

  std::ostringstream help;
  help << "\
Usage: " << argv[0] << " [options] last-train.out alignments.maf\n\
\n\
Find nucleotide substitutions relative to a reference genome.\n\
\n\
Options:\n\
  -h, --help            show this help message and exit\n\
  -b BIAS, --bias=BIAS  require that the strand bias has magnitude < BIAS\n\
  -m INC, --min=INC     minimum increase in log10(likelihood) over homozygous\n\
                        reference (default=" << args.min << ")\n\
  -p N, --ploidy=N      1=haploid, 2=diploid, etc (default='"
       << args.ploidy[0] << "')\n\
  -f BP, --furthest=BP  only use query sequences with colinear alignments\n\
                        separated by <= BP\n\
  -s BP, --splice=BP    only use query sequences with strong splice signals\n\
                        and a splice >= BP\n\
  -S SIZE, --buffer-size=SIZE        memory limit (default="
       << (args.buffer_size / 1024 / 1024 / 1024) << "G)\n\
  -T DIR, --temporary-directory=DIR  put temporary files in DIR\n\
  -v, --verbose         show progress messages\n\
  -V, --version         show version number and exit\n\
";

  const char sOpts[] = "hm:b:p:f:s:S:T:vV";

  static struct option lOpts[] = {
    { "help",                no_argument,       0, 'h' },
    { "min",                 required_argument, 0, 'm' },
    { "bias",                required_argument, 0, 'b' },
    { "ploidy",              required_argument, 0, 'p' },
    { "furthest",            required_argument, 0, 'f' },
    { "splice",              required_argument, 0, 's' },
    { "buffer-size",         required_argument, 0, 'S' },
    { "temporary-directory", required_argument, 0, 'T' },
    { "verbose",             no_argument,       0, 'v' },
    { "version",             no_argument,       0, 'V' },
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
    case 'b':
      args.bias = doubleFromString(optarg);
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
    case 'S':
      args.buffer_size = sizeFromString(optarg);
      break;
    case 'T':
      args.temporary_directory = optarg;
      break;
    case 'v':
      args.verbose++;
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
