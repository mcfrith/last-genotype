#! /bin/sh

# XXX too slow

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=..:$PATH

{
    try last-sub -h
    try last-sub hs-rna.mat hs-rna.maf
    try last-sub -m5 hs-rna.mat hs-rna.maf
    try last-sub -p1 hs-rna.mat hs-rna.maf
    try last-sub -p3 hs-rna.mat hs-rna.maf
} 2>&1 |
diff -u $(basename $0 .sh).out -
