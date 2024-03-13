CXXFLAGS = -O3 -Wall -g

SRC = last-genotype.cc last-genotype-main.cc

INC = last-genotype.hh mcf_string_view.hh mcf_tmpfile.hh	\
mcf_zstream.hh version.hh

last-genotype: ${SRC} ${INC}
	${CXX} ${CXXFLAGS} -o $@ ${SRC} -lz

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e .git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo ${VERSION} | cmp -s $@ - || echo ${VERSION} > $@

FORCE:

clean:
	rm -f last-genotype

tag:
	git tag -m "" `git log --oneline | grep -c .`
