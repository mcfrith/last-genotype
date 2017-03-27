// Copyright 2017 Martin C. Frith

// This routine aims to make a temporary file, in a specified
// directory, and guarantee that the file will be gone when the
// program finishes.

#ifndef MCF_TMPFILE_HH
#define MCF_TMPFILE_HH

#include <cstring>
#include <stdexcept>
#include <stdio.h>  // fdopen
#include <stdlib.h>  // mkstemp
#include <vector>

namespace mcf {

inline FILE *tmpFile(const char *directory) {
  if (!directory) directory = getenv("TMPDIR");
  if (!directory) directory = "/tmp";
  char suffix[] = "/XXXXXX";
  std::vector<char> v(std::strlen(directory) + sizeof suffix);
  char *fileName = &v[0];
  std::strcpy(fileName, directory);
  std::strcat(fileName, suffix);
  int fd = mkstemp(fileName);
  if (fd == -1) throw std::runtime_error("mkstemp temporary file failed");
  FILE *f = fdopen(fd, "wb+");
  if (!f) throw std::runtime_error("fdopen temporary file failed");
  int r = remove(fileName);
  if (r) throw std::runtime_error("remove temporary file failed");
  return f;
}

}

#endif
