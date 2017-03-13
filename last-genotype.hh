// Copyright 2017 Martin C. Frith

#include <vector>

struct LastGenotypeArguments {
  char **argv;
  double min;
  std::vector<const char *> ploidy;
  double furthest;
  double splice;
  const char *lastTrainFile;
  char **mafFiles;
};

void lastGenotype(const LastGenotypeArguments &args);
