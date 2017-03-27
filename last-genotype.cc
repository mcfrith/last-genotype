// Copyright 2017 Martin C. Frith

#include "last-genotype.hh"

#include "mcf_string_view.hh"

#include <fnmatch.h>

#include <cassert>
#include <cstdlib>  // strtod
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>

using std::vector;
using namespace mcf;

typedef unsigned char uchar;
typedef const char *const_char_p;

const size_t maxSize = -1;
const size_t numOfChars = UCHAR_MAX + 1;
const unsigned alphLen = 4;
const unsigned alphLen2 = alphLen * 2;

struct PloidySpec {
  std::string seqNamePattern;
  unsigned ploidy;
};

struct Alignment {
  unsigned refSeqNum;
  unsigned beg;
  unsigned end;
  unsigned querySeqNum;  // xxx size_t ?
  uchar *columns;
};

struct AlignedBase {
  unsigned querySeqNum;
  uchar queryBase;
  double prob;
};

struct AlignedBaseText {
  char t[4];
};

static unsigned alignmentStrandNum(const Alignment &a) {
  return a.columns[0] % 2;
}

static const uchar *columnFromAlignment(const Alignment &a, size_t coord) {
  return a.columns + (coord - a.beg) * 3;
}

static size_t alignmentDistance(const Alignment &a, const Alignment &b) {
  return a.end <= b.beg ? b.beg - a.end : maxSize;
}

static bool isLessByGenome(const Alignment &a, const Alignment &b) {
  if (a.refSeqNum != b.refSeqNum) return a.refSeqNum < b.refSeqNum;
  if (a.beg != b.beg) return a.beg < b.beg;
  // xxx stuff below here unnecessary?
  if (a.end != b.end) return a.end < b.end;
  return std::memcmp(a.columns, b.columns, a.end - a.beg) < 0;
}

static bool isLessByQuery(const AlignedBase &a, const AlignedBase &b) {
  return a.querySeqNum < b.querySeqNum;
}

static bool isLessByAscii(const AlignedBaseText &a, const AlignedBaseText &b) {
  return std::strcmp(a.t, b.t) < 0;
}

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static StringView &munch(StringView &in, const char *text) {
  // xxx skip whitespace?
  const char *b = in.begin();
  const char *e = in.end();
  while (*text) {
    if (b == e || *b != *text) return in = StringView();
    ++text;
    ++b;
  }
  return in = StringView(b, e);
}

static std::istream &openIn(const char *fileName, std::ifstream &ifs) {
  if (isChar(fileName, '-')) return std::cin;
  ifs.open(fileName);
  if (!ifs) err("can't open file: " + std::string(fileName));
  return ifs;
}

static double myLog(double x) {
  return log(x);  // change to e.g. log2 if faster
}

static void makeQualTable(double *qualTable) {  // fastq-sanger ascii codes
  size_t offset = 33;
  for (size_t i = 0; i < numOfChars; ++i)
    qualTable[i] = (i < offset) ? 1 : 1 - pow(10.0, -0.1 * (i - offset));
}

static void ploidiesFromStrings(const const_char_p *argsBeg,
				const const_char_p *argsEnd,
				vector<PloidySpec> &ploidies) {
  for (const const_char_p *i = argsBeg; i < argsEnd; ++i) {
    for (const char *beg = *i; ; ) {
      const char *end = beg + std::strcspn(beg, ",");
      const char *mid = std::find(beg, end, ':');
      PloidySpec spec;
      if (mid < end) {
	spec.seqNamePattern.assign(beg, mid);
	beg = mid + 1;
      } else {
	spec.seqNamePattern = "*";
      }
      StringView s(beg, end);
      s >> spec.ploidy;
      if(!s || !s.empty() || spec.ploidy < 1) err("bad ploidy");
      ploidies.push_back(spec);
      if (*end == 0) break;
      beg = end + 1;
    }
  }
  reverse(ploidies.begin(), ploidies.end());
}

static unsigned ploidyOfChromosome(const vector<PloidySpec> &ploidies,
				   const char *chromosomeName) {
  for (size_t i = 0; i < ploidies.size(); ++i)
    // xxx FNM_NOESCAPE?
    if (fnmatch(ploidies[i].seqNamePattern.c_str(), chromosomeName, 0) == 0)
      return ploidies[i].ploidy;
  assert(0);
}

// the number of genotypes with alphabet size a and ploidy p:
// = (a + p - 1) choose p
// = (a + p - 1)! / [p! (a - 1)!]
// = (p+1) * (p+2) * (p+3) / 6

static void makeGenotypes(vector<uchar> &genotypes, unsigned ploidy) {
  genotypes.assign(ploidy, 0);  // start with the all-A genotype
  while (true) {
    size_t end = genotypes.size();  // end of previous genotype
    size_t beg = end - ploidy;      // beginning of previous genotype
    size_t mid = end;
    unsigned x = alphLen;
    while (mid > beg) {
      --mid;
      x = genotypes[mid] + 1;
      if (x < alphLen) break;
    }
    if (x == alphLen) break;
    for (; beg < mid; ++beg) genotypes.push_back(genotypes[beg]);
    for (; mid < end; ++mid) genotypes.push_back(x);
  }
}

static size_t homozygousGenotypeIndex(unsigned ploidy, unsigned baseCode) {
  size_t r = 0;
  switch (baseCode) {
  case 3:
    r += ploidy;
  case 2:
    r += (ploidy + 1) * ploidy / 2;
  case 1:
    r += (ploidy + 2) * (ploidy + 1) * ploidy / 6;
  }
  return r;
}

static void makeGenotypeCalc(double baseCalcMatrix[][alphLen2],
			     unsigned ploidy,
			     const vector<uchar> &genotypes,
			     vector<double> &genotypeCalcMatrix) {
  size_t numOfGenotypes = genotypes.size() / ploidy;
  genotypeCalcMatrix.resize(alphLen2 * numOfGenotypes);
  for (size_t i = 0; i < numOfGenotypes; ++i) {
    for (unsigned j = 0; j < alphLen2; ++j) {
      double sum = 0;
      for (unsigned k = 0; k < ploidy; ++k) {
	unsigned genomicBase = genotypes[ploidy * i + k];
	sum += baseCalcMatrix[genomicBase][j];
      }
      genotypeCalcMatrix[alphLen2 * i + j] = sum / ploidy;
    }
  }
}

static void alignmentProbMatrixFromLastTrain(const char *fileName,
					     double matrix[][alphLen]) {
  std::ifstream ifs;
  std::istream &in = openIn(fileName, ifs);

  unsigned rowNum = alphLen + 1;
  for (std::string line; getline(in, line); ) {
    if (rowNum >= alphLen) {
      StringView s(line);
      if (munch(s, "# probability matrix")) rowNum = 0;
    } else {
      std::istringstream s(line);
      std::string junk;
      s >> junk >> junk;
      for (unsigned i = 0; i < alphLen; ++i) s >> matrix[rowNum][i];
      if (s) ++rowNum;
    }
  }
  if (rowNum != alphLen) err("can't read last-train file");
}

static void baseCalcMatrixFromLastTrain(const char *fileName,
					double baseCalcMatrix[][alphLen2]) {
  double matrix[alphLen][alphLen];
  alignmentProbMatrixFromLastTrain(fileName, matrix);

  // Here, we replace prob(i and j) with: prob(j given i) - 1.
  // The "- 1" is cryptic, but it simplifies some later calculations.
  for (unsigned i = 0; i < alphLen; ++i) {
    double sum = 0;
    for (unsigned j = 0; j < alphLen; ++j) {
      if (matrix[i][j] < 0) err("negative last-train probability");
      sum += matrix[i][j];
    }
    if (sum <= 0) err("bad last-train probabilities");
    for (unsigned j = 0; j < alphLen; ++j) {
      matrix[i][j] = matrix[i][j] / sum - 1;
    }
  }

  for (unsigned i = 0; i < alphLen; ++i) {
    for (unsigned j = 0; j < alphLen; ++j) {
      baseCalcMatrix[i][j * 2] = matrix[i][j];
      baseCalcMatrix[i][j * 2 + 1] = matrix[alphLen - i - 1][alphLen - j - 1];
    }
  }
}

static void calcMinCoverage(double minLogProbIncrease,
			    double baseCalcMatrix[][alphLen2],
			    double *minCoveragePerRefBase) {
  for (unsigned r = 0; r < alphLen; ++r) {  // loop over reference ACGT
    double maxRatio = 1;
    for (unsigned q = 0; q < alphLen2; ++q) {
      double maxProb = 0;
      for (unsigned g = 0; g < alphLen; ++g) {
	double probQueryGivenRef = baseCalcMatrix[g][q] + 1;
	maxProb = std::max(maxProb, probQueryGivenRef);
      }
      if (maxProb <= 0) continue;
      double probQueryGivenRef = baseCalcMatrix[r][q] + 1;
      double ratio = maxProb / probQueryGivenRef;  // xxx may be div 0 -> inf
      maxRatio = std::max(maxRatio, ratio);
    }
    minCoveragePerRefBase[r] = minLogProbIncrease / myLog(maxRatio);
    // xxx may be div 0 -> inf
  }

  std::cout << "# Minimum coverage per reference base:";
  for (unsigned i = 0; i < alphLen; ++i) {
    std::cout << ' ' << "ACGT"[i] << '=' << ceil(minCoveragePerRefBase[i]);
  }
  std::cout << '\n';
}

static size_t stringIndex(vector<std::string> &strings, StringView s) {
  // For now, this uses trivial linear search.  If there are many
  // strings, binary search might be better.
  size_t n = strings.size();
  for (size_t i = 0; i < n; ++i) {
    StringView t(strings[i]);
    if (t == s) return i;
  }
  strings.push_back(std::string(s.begin(), s.size()));
  return n;
}

static void makeSeqCodeTables(uchar seqCodeTables[][numOfChars]) {
  for (size_t i = 0; i < numOfChars; ++i) {
    seqCodeTables[0][i] = alphLen * 16;
    seqCodeTables[1][i] = alphLen * 2;
    seqCodeTables[2][i] = alphLen * 2 + 1;
  }

  for (unsigned i = 0; i < alphLen; ++i) {
    uchar x = "ACGT"[i];
    uchar y = "acgt"[i];
    seqCodeTables[0][x] = seqCodeTables[0][y] = i * 16;
    seqCodeTables[1][x] = seqCodeTables[1][y] = i * 2;
    seqCodeTables[2][x] = seqCodeTables[2][y] = i * 2 + 1;
  }
}

static size_t basesBetweenAlignments(const Alignment &a, const Alignment &b) {
  if (a.refSeqNum != b.refSeqNum) return maxSize;
  unsigned s = alignmentStrandNum(a);
  if (s != alignmentStrandNum(b)) return maxSize;
  return (s == 0) ? alignmentDistance(a, b) : alignmentDistance(b, a);
}

static bool isGoodAlignments(const LastGenotypeArguments &args,
			     const Alignment *beg, const Alignment *end) {
  assert(beg < end);
  bool isGood = (args.splice < 0);
  for ( ; beg < end - 1; ++beg) {
    size_t d = basesBetweenAlignments(beg[0], beg[1]);
    if (args.furthest >= 0) {
      if (d == maxSize || d > args.furthest) return false;
    }
    if (args.splice >= 0) {
      if (d < maxSize && d >= args.splice) isGood = true;
    }
  }
  return isGood;
}

static void parseHeadLine(const std::string &aLine, double &sense) {
  StringView s(aLine);
  StringView w;
  while (s >> w) {
    if (munch(w, "sense=")) {
      // xxx ugly
      std::string value(w.data(), w.size());
      const char *b = value.c_str();
      char *e;
      sense = std::strtod(b, &e);
      // don't bother with errno
      if (e > b) return;
    }
  }
  err("missing sense=");
}

static void parseTopSeq(const std::string &sLine, StringView &seqName,
			unsigned &start, StringView &seq) {
  StringView s(sLine);
  s >> seq >> seqName >> start >> seq >> seq >> seq >> seq;
  if (!s) err("bad MAF line: " + sLine);
}

static void parseBotSeq(const std::string &sLine, StringView &seqName,
			StringView &strand, StringView &seq) {
  StringView s(sLine);
  s >> seq >> seqName >> seq >> seq >> strand >> seq >> seq;
  if (!s) err("bad MAF line: " + sLine);
}

static void parseProbLine(const std::string &pLine, StringView &probSeq) {
  StringView s(pLine);
  s >> probSeq >> probSeq;
  if (!s) err("bad MAF line: " + pLine);
}

static bool isBadQuery(const LastGenotypeArguments &args,
		       const std::string &aLine) {
  if (args.splice >= 0) {
    double sense;
    parseHeadLine(aLine, sense);
    return fabs(sense) < 10;  // xxx hard-coded empirical value
  }
  return false;
}

static void doOneQuery(const LastGenotypeArguments &args,
		       vector<Alignment> &alignments,
		       size_t &queryStart, size_t &queryCount) {
  size_t s = alignments.size();
  if (queryStart >= s) return;
  const Alignment *a = &alignments[0];
  if (isGoodAlignments(args, a + queryStart, a + s)) {
    queryStart = s;
    ++queryCount;
  } else {
    while (s > queryStart) {
      --s;
      delete[] alignments[s].columns;
    }
    alignments.resize(s);
  }
}

static size_t alignmentSpan(StringView sequence) {
  return sequence.size() - std::count(sequence.begin(), sequence.end(), '-');
}

static void checkSequenceLength(StringView sequence, size_t length) {
  if (sequence.size() != length) err("unequal MAF block lengths");
}

static uchar *alignmentColumns(uchar seqCodeTables[][numOfChars],
			       size_t colBytes,
			       StringView strand,
			       StringView rSeq,
			       StringView qSeq,
			       const StringView *probSeqs,
			       unsigned pLineCount) {
  size_t n = rSeq.size();
  checkSequenceLength(qSeq, n);
  if (pLineCount > 0) checkSequenceLength(probSeqs[0], n);
  if (pLineCount > 1) checkSequenceLength(probSeqs[1], n);
  const uchar *rTable = seqCodeTables[0];
  const uchar *qTable = seqCodeTables[(strand[0] == '-') ? 2 : 1];
  uchar *columns = new uchar[colBytes];
  uchar *c = columns;
  for (size_t i = 0; i < n; ++i) {
    uchar r = rSeq[i];
    if (r == '-') continue;
    uchar q = qSeq[i];
    *c++ = rTable[r] + qTable[q];
    *c++ = (pLineCount > 0) ? probSeqs[0][i] : 0;
    *c++ = (pLineCount > 1) ? probSeqs[1][i] : 0;
  }
  return columns;
}

static void addAlignment(vector<Alignment> &alignments,
			 size_t refSeqNum, unsigned beg, unsigned len,
			 size_t querySeqNum, uchar *columns) {
  Alignment a;
  a.refSeqNum = refSeqNum;
  a.beg = beg;
  a.end = beg + len;
  a.querySeqNum = querySeqNum;
  a.columns = columns;
  if (a.refSeqNum < refSeqNum) err("too many reference sequences");
  if (a.querySeqNum < querySeqNum) err("too many query sequences");
  alignments.push_back(a);
}

static void readMaf(const LastGenotypeArguments &args,
		    vector<Alignment> &alignments,
		    vector<std::string> &refSeqNames,
		    size_t &queryCount,
		    std::istream &in) {
  uchar seqCodeTables[3][numOfChars];
  makeSeqCodeTables(seqCodeTables);

  bool isBad = false;
  size_t queryStart = alignments.size();
  unsigned sLineCount = 0;
  unsigned pLineCount = 0;
  std::string aLine;
  std::string sLineBuf[4];
  std::string *sLines = sLineBuf;
  std::string pLines[2];
  std::string line;
  StringView rName, qName, qNameOld, strand, rSeq, qSeq, probSeqs[2];

  do {
    getline(in, line);
    const char *s = line.c_str();
    if (*s == 'a') {
      line.swap(aLine);
    } else if (*s == 's') {
      if (sLineCount > 1) err("too many MAF s lines");
      line.swap(sLines[sLineCount]);
      ++sLineCount;
    } else if (*s == 'p') {
      if (pLineCount > 1) err("too many MAF p lines");
      line.swap(pLines[pLineCount]);
      parseProbLine(pLines[pLineCount], probSeqs[pLineCount]);
      ++pLineCount;
    } else if (!isGraph(*s)) {
      if (sLineCount > 1) {
	unsigned rBeg;
	parseTopSeq(sLines[0], rName, rBeg, rSeq);
	parseBotSeq(sLines[1], qName, strand, qSeq);
	if (qName != qNameOld) {
	  doOneQuery(args, alignments, queryStart, queryCount);
	  isBad = isBadQuery(args, aLine);
	}
	if (!isBad) {
	  size_t refSeqNum = stringIndex(refSeqNames, rName);
	  size_t rLen = alignmentSpan(rSeq);
	  size_t colBytes = rLen * 3;
	  uchar *columns = alignmentColumns(seqCodeTables, colBytes, strand,
					    rSeq, qSeq, probSeqs, pLineCount);
	  addAlignment(alignments, refSeqNum, rBeg, rLen, queryCount, columns);
	}
	sLines = sLineBuf + (sLines - sLineBuf + 2) % 4;
	qNameOld = qName;
      }
      sLineCount = pLineCount = 0;
    }
  } while (in);

  doOneQuery(args, alignments, queryStart, queryCount);
}

static void readAlignmentFiles(const LastGenotypeArguments &args,
			       vector<Alignment> &alignments,
			       vector<std::string> &refSeqNames) {
  size_t queryCount = 0;

  if (*args.mafFiles) {
    for (char **i = args.mafFiles; *i; ++i) {
      std::ifstream ifs;
      std::istream &in = openIn(*i, ifs);
      readMaf(args, alignments, refSeqNames, queryCount, in);
    }
  } else {
    readMaf(args, alignments, refSeqNames, queryCount, std::cin);
  }

  sort(alignments.begin(), alignments.end(), isLessByGenome);

  std::cout << "# Query sequences used: " << queryCount << '\n';
  std::cout << "# Alignments used: " << alignments.size() << '\n';
}

static size_t preprocessColumns(const double *qualTable, size_t coord,
				const vector<Alignment> &alignments,
				vector<uchar> &colBases,
				vector<double> &colProbs) {
  size_t n = alignments.size();
  if (colBases.size() < n) colBases.resize(n);
  if (colProbs.size() < n) colProbs.resize(n);
  size_t numOfBases = 0;
  for (size_t i = 0; i < n; ++i) {
    const uchar *c = columnFromAlignment(alignments[i], coord);
    uchar queryBase = c[0] % 16;
    if (queryBase >= alphLen2) continue;
    colBases[numOfBases] = queryBase;
    colProbs[numOfBases] = qualTable[c[1]] * qualTable[c[2]];
    ++numOfBases;
  }
  return numOfBases;
}

static void getAlignedBases(const double *qualTable, size_t coord,
			    const vector<Alignment> &alignments,
			    vector<AlignedBase> &alignedBases) {
  size_t j = 0;
  for (size_t i = 0; i < alignments.size(); ++i) {
    const uchar *c = columnFromAlignment(alignments[i], coord);
    uchar queryBase = c[0] % 16;
    if (queryBase >= alphLen2) continue;
    alignedBases[j].querySeqNum = alignments[i].querySeqNum;
    alignedBases[j].queryBase = queryBase;
    alignedBases[j].prob = qualTable[c[1]] * qualTable[c[2]];
    ++j;
  }
}

static void makeAlignedBaseTexts(size_t coord,
				 const vector<Alignment> &alignments,
				 vector<AlignedBaseText> &alignedBaseTexts) {
  size_t j = 0;
  for (size_t i = 0; i < alignments.size(); ++i) {
    const uchar *c = columnFromAlignment(alignments[i], coord);
    uchar queryBase = c[0] % 16;
    if (queryBase >= alphLen2) continue;
    alignedBaseTexts[j].t[0] = "AaCcGgTt"[queryBase];
    alignedBaseTexts[j].t[1] = c[1];
    alignedBaseTexts[j].t[2] = c[2];
    ++j;
  }
}

static bool isAll(unsigned refBase, size_t numOfBases, const uchar *colBases) {
  for (size_t i = 0; i < numOfBases; ++i)
    if (colBases[i] / 2 != refBase) return false;
  return true;
}

static void calcLogProbs(const double *genotypeCalcMatrix,
			 size_t numOfBases,
			 const uchar *colBases,
			 const double *colProbs,
			 vector<double> &genotypeLogProbs) {
  size_t numOfGenotypes = genotypeLogProbs.size();
  for (size_t i = 0; i < numOfGenotypes; ++i) {
    const double *g = &genotypeCalcMatrix[alphLen2 * i];
    double logProb = 0;
    for (size_t j = 0; j < numOfBases; ) {
      // Tests on chimera h007:
      // last-genotype -f1e6 -s50 lc2adbulk2-pass_2d.mat
      //                          lc2adbulk2-pass_2d-d90-m50-j4.maf
      // Times without batching:
      // log10: about 12 seconds
      // log: about 11 seconds
      // log2: about 7.6 seconds
      // log2f: about 6.6 seconds
      // log1p: about 6.2 seconds
      // log1pf: about 6.3 seconds
      // With batching: got overflow with batchSize 256, but not 128
      // batchSize 64 was reproducibly faster than 32
      const unsigned batchSize = 64;
      size_t end = std::min(j + batchSize, numOfBases);
      double prob = 1;
      for (; j < end; ++j) {
	prob *= colProbs[j] * g[colBases[j]] + 1;
      }
      logProb += myLog(prob);  // xxx could be log 0, which should be -inf
    }
    genotypeLogProbs[i] = logProb;
  }
}

static double twoSiteLogProb(size_t numOfPairedBases,
			     const uchar *pairBases, const double *pairProbs,
			     const double *rowA1, const double *rowA2,
			     const double *rowB1, const double *rowB2) {
  // 1, 2 refer to the two sites
  // A, B refer to the maternal & paternal chromosomes
  double logProb = 0;
  for (size_t i = 0; i < numOfPairedBases; ++i) {
    uchar base1 = pairBases[i];  // the query/read base at site 1
    double prob1 = pairProbs[i];  // probability that it's correctly aligned
    ++i;
    uchar base2 = pairBases[i];  // the query/read base at site 2
    double prob2 = pairProbs[i];  // probability that it's correctly aligned
    double a = (prob1 * rowA1[base1] + 1) * (prob2 * rowA2[base2] + 1);
    double b = (prob1 * rowB1[base1] + 1) * (prob2 * rowB2[base2] + 1);
    logProb += myLog(0.5 * (a + b));  // xxx we could have log(0) -> -inf
  }
  return logProb;
}

static void findTopTwo(const vector<double> &v, size_t &max1, size_t &max2) {
  size_t s = v.size();
  assert(s > 1);
  max1 = 0;
  for (size_t i = 1; i < s; ++i) {
    if (v[i] > v[max1]) max1 = i;
  }
  max2 = max1 ? 0 : 1;
  for (size_t i = 1; i < s; ++i) {
    if (v[i] > v[max2] && i != max1) max2 = i;
  }
}

static size_t findPairs(const vector<AlignedBase> &x,
			const vector<AlignedBase> &y,
			vector<uchar> &pairBases, vector<double> &pairProbs) {
  size_t xs = x.size();
  size_t ys = y.size();
  if (pairBases.size() < xs + ys) pairBases.resize(xs + ys);
  if (pairProbs.size() < xs + ys) pairProbs.resize(xs + ys);
  size_t xi = 0;
  size_t yi = 0;
  size_t j = 0;
  while (xi < xs && yi < ys) {
    size_t xq = x[xi].querySeqNum;
    size_t yq = y[yi].querySeqNum;
    if (xq < yq) {
      ++xi;
    } else if (yq < xq) {
      ++yi;
    } else if ((xi + 1 < xs && x[xi + 1].querySeqNum == xq) ||
	       (yi + 1 < ys && y[yi + 1].querySeqNum == yq)) {
      do { ++xi; } while (xi < xs && x[xi].querySeqNum == xq);
      do { ++yi; } while (yi < ys && y[yi].querySeqNum == yq);
    } else {
      pairBases[j] = x[xi].queryBase;
      pairProbs[j] = x[xi].prob;
      ++xi;
      ++j;
      pairBases[j] = y[yi].queryBase;
      pairProbs[j] = y[yi].prob;
      ++yi;
      ++j;
    }
  }
  return j;
}

static double doPhase(double baseCalcMatrix[][alphLen2],
		      const vector<uchar> &genotype1, vector<uchar> &genotype2,
		      size_t numOfPairedBases,
		      const uchar *pairBases, const double *pairProbs) {
  // 1, 2 refer to the two sites
  // A, B refer to the maternal & paternal chromosomes
  if (genotype1.size() != 2) return 0;
  const double *rowA1 = baseCalcMatrix[genotype1[0]];
  const double *rowA2 = baseCalcMatrix[genotype2[0]];
  const double *rowB1 = baseCalcMatrix[genotype1[1]];
  const double *rowB2 = baseCalcMatrix[genotype2[1]];
  double p = twoSiteLogProb(numOfPairedBases, pairBases, pairProbs,
			    rowA1, rowA2, rowB1, rowB2);
  double q = twoSiteLogProb(numOfPairedBases, pairBases, pairProbs,
			    rowA1, rowB2, rowB1, rowA2);
  if (p >= q) return p - q;
  std::reverse(genotype2.begin(), genotype2.end());
  return q - p;
}

static void decodeGenotype(unsigned ploidy, const uchar *in, char *out) {
  for (unsigned i = 0; i < ploidy; ++i) {
    assert(in[i] < alphLen);
    out[i] = "ACGT"[in[i]];
  }
}

void discardOldAlignments(vector<Alignment> &alignments, size_t coord) {
  size_t n = alignments.size();
  size_t j = 0;
  for (size_t i = 0; i < n; ++i) {
    if (alignments[i].end > coord) {
      alignments[j++] = alignments[i];
    }
  }
  alignments.resize(j);
}

static void printArgs(char **argv) {
  std::cout << "# last-genotype version "
#include "version.hh"
    "\n";
  std::cout << '#';
  for (char **i = argv; *i; ++i) std::cout << ' ' << *i;
  std::cout << '\n';
}

void lastGenotype(const LastGenotypeArguments &args) {
  double qualTable[numOfChars];
  makeQualTable(qualTable);

  vector<PloidySpec> ploidies;
  ploidiesFromStrings(&args.ploidy.front(), &args.ploidy.back() + 1, ploidies);

  printArgs(args.argv);

  double baseCalcMatrix[alphLen][alphLen2];
  baseCalcMatrixFromLastTrain(args.lastTrainFile, baseCalcMatrix);

  double minLogProbIncrease = args.min * myLog(10);
  double minCoveragePerRefBase[alphLen];
  calcMinCoverage(minLogProbIncrease, baseCalcMatrix, minCoveragePerRefBase);

  vector<Alignment> alignments;
  vector<std::string> refSeqNames;
  readAlignmentFiles(args, alignments, refSeqNames);

  size_t numOfTestedSites = 0;
  size_t alignmentPos = 0;
  size_t refSeqNum = 0;
  size_t coord = 0;
  unsigned ploidy = 0;  // shut the compiler up
  vector<uchar> genotypes;
  vector<double> genotypeCalcMatrix;
  vector<Alignment> alignmentsHere;
  vector<double> genotypeLogProbs;
  vector<uchar> colBases;
  vector<double> colProbs;
  vector<AlignedBase> oldAlignedBases, newAlignedBases;
  vector<uchar> oldGenotype, newGenotype;
  vector<char> genotypeString;
  vector<AlignedBaseText> alignedBaseTexts;

  std::cout.precision(2);

  while (true) {
    if (alignmentsHere.empty()) {
      if (alignmentPos == alignments.size()) break;
      const Alignment &a = alignments[alignmentPos];
      refSeqNum = a.refSeqNum;
      coord = a.beg;
      ploidy = ploidyOfChromosome(ploidies, refSeqNames[refSeqNum].c_str());
      makeGenotypes(genotypes, ploidy);
      makeGenotypeCalc(baseCalcMatrix, ploidy, genotypes, genotypeCalcMatrix);
      genotypeLogProbs.resize(genotypes.size() / ploidy);
      newGenotype.resize(ploidy);
      genotypeString.assign(ploidy + 1, 0);
    }
    while (alignmentPos < alignments.size()) {
      const Alignment &a = alignments[alignmentPos];
      if (a.refSeqNum > refSeqNum || a.beg > coord) break;
      alignmentsHere.push_back(a);
      ++alignmentPos;
    }
    while (true) {  // xxx bogus loop: what's the right way to do this?
      unsigned refBase = columnFromAlignment(alignmentsHere[0], coord)[0] / 16;
      if (refBase >= alphLen) break;
      if (alignmentsHere.size() < minCoveragePerRefBase[refBase]) break;
      ++numOfTestedSites;
      size_t numOfBases = preprocessColumns(qualTable, coord, alignmentsHere,
					    colBases, colProbs);
      if (isAll(refBase, numOfBases, &colBases[0])) break;  // makes it faster
      calcLogProbs(&genotypeCalcMatrix[0], numOfBases,
		   &colBases[0], &colProbs[0], genotypeLogProbs);
      size_t refGenotypeIndex = homozygousGenotypeIndex(ploidy, refBase);
      double refLogProb = genotypeLogProbs[refGenotypeIndex];
      size_t max1, max2;
      findTopTwo(genotypeLogProbs, max1, max2);
      double logProb1 = genotypeLogProbs[max1];
      double logProbIncRef = logProb1 - refLogProb;
      if (logProbIncRef < minLogProbIncrease) break;
      double logProb2 = genotypeLogProbs[max2];
      double logProbInc2nd = logProb1 - logProb2;
      std::memcpy(&newGenotype[0], &genotypes[max1 * ploidy], ploidy);
      const uchar *genotype2nd = &genotypes[max2 * ploidy];
      double logProbIncPhase = 0;
      size_t phaseCoverage = 0;
      if (ploidy == 2 && newGenotype[0] != newGenotype[1]) {
	newAlignedBases.resize(numOfBases);
	getAlignedBases(qualTable, coord, alignmentsHere, newAlignedBases);
	sort(newAlignedBases.begin(), newAlignedBases.end(), isLessByQuery);
	size_t numOfPairedBases = findPairs(oldAlignedBases, newAlignedBases,
					    colBases, colProbs);
	logProbIncPhase = doPhase(baseCalcMatrix, oldGenotype, newGenotype,
				  numOfPairedBases,
				  &colBases[0], &colProbs[0]);
	phaseCoverage = numOfPairedBases / 2;
	oldAlignedBases.swap(newAlignedBases);
	oldGenotype = newGenotype;
      }

      alignedBaseTexts.resize(numOfBases);
      makeAlignedBaseTexts(coord, alignmentsHere, alignedBaseTexts);
      sort(alignedBaseTexts.begin(), alignedBaseTexts.end(), isLessByAscii);

      std::cout << refSeqNames[refSeqNum] << '\t'
		<< coord << '\t'
		<< "ACGT"[refBase] << '\t';
      decodeGenotype(ploidy, &newGenotype[0], &genotypeString[0]);
      std::cout << &genotypeString[0] << '\t'
		<< (logProbIncRef / myLog(10)) << '\t'
		<< (logProbInc2nd / myLog(10)) << '\t';
      decodeGenotype(ploidy, genotype2nd, &genotypeString[0]);
      std::cout << &genotypeString[0] << '\t'
		<< (logProbIncPhase / myLog(10)) << '\t'
		<< phaseCoverage << '\t';
      for (size_t i = 0; i < numOfBases; ++i) {
	if (i) std::cout << ' ';
	std::cout << alignedBaseTexts[i].t;
      }
      std::cout << '\n';
      break;
    }
    ++coord;
    discardOldAlignments(alignmentsHere, coord);
  }

  std::cout << "# Tested sites: " << numOfTestedSites << '\n';
}
