// Copyright 2016 Martin C. Frith
// Yet another re-invention of basic text parsing functions.
// Incomplete: add more functions when needed.

#ifndef MCF_STRING_VIEW_HH
#define MCF_STRING_VIEW_HH

#include <algorithm>
#include <climits>
#include <cstring>
#include <ostream>
#include <string>
#include <stddef.h>

namespace mcf {

inline bool isGraph(char c) {
  return c > ' ';  // faster than std::isgraph
}

inline bool isDigit(char c) {
  return c >= '0' && c <= '9';
}

inline bool isChar(const char *myString, char myChar) {
  return myString[0] == myChar && myString[1] == 0;
}

class StringView {
public:
  StringView() : b(0), e(0) {}

  explicit StringView(const std::string &s)
    : b(s.data()), e(s.data() + s.size()) {}

  StringView(const char *beg, const char *end) : b(beg), e(end) {}

  operator const void *() const { return b; }

  const char *data()  const { return b; }
  const char *begin() const { return b; }
  const char *end()   const { return e; }
  size_t size()       const { return e - b; }
  bool empty()        const { return b == e; }

  char operator[](size_t i) const { return b[i]; }

  char front() const { return begin()[0]; }
  char back()  const { return end()[-1]; }

  void remove_prefix(size_t n) { b += n; }
  void remove_suffix(size_t n) { e -= n; }

  int compare(StringView v) const {
    size_t len = size();
    size_t vlen = v.size();
    int c = std::memcmp(data(), v.data(), std::min(len, vlen));
    return c ? c : len == vlen ? 0 : len < vlen ? -1 : 1;
  }

private:
  const char *b;
  const char *e;
};

inline bool isDigit(StringView s) {
  const char *b = s.begin();
  const char *e = s.end();
  if (b == e) return false;
  while (1) {
    if (!isDigit(*b)) return false;
    ++b;
    if (b == e) return true;
  }
}

inline bool operator==(StringView x, StringView y) {
  return
    x.size() == y.size() && std::memcmp(x.data(), y.data(), x.size()) == 0;
}

inline bool operator!=(StringView x, StringView y) {
  return !(x == y);
}

inline bool operator< (StringView x, StringView y) {
  return x.compare(y) <  0;
}

inline bool operator<=(StringView x, StringView y) {
  return x.compare(y) <= 0;
}

inline bool operator> (StringView x, StringView y) {
  return x.compare(y) >  0;
}

inline bool operator>=(StringView x, StringView y) {
  return x.compare(y) >= 0;
}

inline bool operator==(StringView x, char y) {
  return x.size() == 1 && x[0] == y;
}

inline bool operator!=(StringView x, char y) {
  return !(x == y);
}

inline bool operator==(StringView x, const char *y) {
  size_t xlen = x.size();
  return std::strncmp(x.data(), y, xlen) == 0 && std::strlen(y) == xlen;
}

inline bool operator!=(StringView x, const char *y) {
  return !(x == y);
}

inline std::ostream &operator<<(std::ostream &out, StringView s) {
  return out.write(s.data(), s.size());
}

inline StringView &operator>>(StringView &in, char &out) {
  const char *b = in.begin();
  const char *e = in.end();
  while (true) {
    if (b == e) return in = StringView();
    if (isGraph(*b)) break;
    ++b;
  }
  out = *b++;
  return in = StringView(b, e);
}

inline StringView &operator>>(StringView &in, StringView &out) {
  const char *b = in.begin();
  const char *e = in.end();
  while (true) {
    if (b == e) return in = StringView();
    if (isGraph(*b)) break;
    ++b;
  }
  const char *m = b;
  do { ++m; } while (m < e && isGraph(*m));
  out = StringView(b, m);
  return in = StringView(m, e);
}

inline StringView &operator>>(StringView &in, unsigned &out) {
  const char *b = in.begin();
  const char *e = in.end();
  while (true) {
    if (b == e) return in = StringView();
    if (isGraph(*b)) break;
    ++b;
  }
  // should we allow an initial '+'?
  if (!isDigit(*b)) return in = StringView();
  unsigned z = *b++ - '0';
  while (b < e && isDigit(*b)) {
    if (z > UINT_MAX / 10) return in = StringView();
    z *= 10;
    unsigned digit = *b++ - '0';
    if (z > UINT_MAX - digit) return in = StringView();
    z += digit;
  }
  out = z;
  return in = StringView(b, e);
}

inline StringView &operator>>(StringView &in, long &out) {
  const char *b = in.begin();
  const char *e = in.end();
  while (true) {
    if (b == e) return in = StringView();
    if (isGraph(*b)) break;
    ++b;
  }
  if (*b == '-') {
    ++b;
    if (b == e || !isDigit(*b)) return in = StringView();
    long z = '0' - *b++;
    while (b < e && isDigit(*b)) {
      if (z < LONG_MIN / 10) return in = StringView();
      z *= 10;
      long digit = *b++ - '0';
      if (z < LONG_MIN + digit) return in = StringView();
      z -= digit;
    }
    out = z;
  } else {
    // should we allow an initial '+'?
    if (!isDigit(*b)) return in = StringView();
    long z = *b++ - '0';
    while (b < e && isDigit(*b)) {
      if (z > LONG_MAX / 10) return in = StringView();
      z *= 10;
      long digit = *b++ - '0';
      if (z > LONG_MAX - digit) return in = StringView();
      z += digit;
    }
    out = z;
  }
  return in = StringView(b, e);
}

}

#endif
