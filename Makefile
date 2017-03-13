CXX = g++
CXXFLAGS = -O3 -Wall -g

SRC = last-genotype.cc last-genotype-main.cc

last-genotype: ${SRC} last-genotype.hh version.hh
	${CXX} ${CXXFLAGS} -o $@ ${SRC}

VERSION = \"`git log --oneline | grep -c .``git diff --quiet HEAD || echo +`\"
UNKNOWN = \"UNKNOWN\"

version.hh: FORCE
	if test -e .git ; \
	then echo ${VERSION} | cmp -s $@ - || echo ${VERSION} > $@ ; \
	else test -e $@ || echo ${UNKNOWN} > $@ ; \
	fi

FORCE:

clean:
	rm -f last-genotype
