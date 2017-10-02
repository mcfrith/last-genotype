CXXFLAGS = -O3 -Wall -g

SRC = last-genotype.cc last-genotype-main.cc

INC = last-genotype.hh mcf_string_view.hh mcf_tmpfile.hh	\
mcf_zstream.hh version.hh

last-genotype: ${SRC} ${INC}
	${CXX} ${CXXFLAGS} -o $@ ${SRC} -lz

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
