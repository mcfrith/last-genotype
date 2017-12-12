// Copyright 2017 Martin C. Frith

#include <stddef.h>
#include <vector>

struct LastGenotypeArguments {
  char **argv;
  double min_ref;
  double min_2nd;
  double bias_ref;
  double bias_2nd;
  std::vector<const char *> ploidy;
  double furthest;
  double splice;
  const char *class_file;
  size_t buffer_size;
  const char *temporary_directory;
  int verbose;
  const char *lastTrainFile;
  char **mafFiles;
};

void lastGenotype(const LastGenotypeArguments &args);
