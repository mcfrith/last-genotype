#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=..:$PATH

{
    try last-genotype -h
    try last-genotype hs-rna.mat hs-rna.maf
    try last-genotype -m5 hs-rna.mat hs-rna.maf
    try last-genotype -p1 hs-rna.mat hs-rna.maf
    try last-genotype -p3 hs-rna.mat hs-rna.maf
    try last-genotype -f1e6 -s50 hs-rna.mat hs-rna.maf
    try last-genotype -f1000 hs-rna.mat hs-rna.maf
    try last-genotype -p'3,chr1*:1' hs-rna.mat hs-rna.maf
    try last-genotype -f1000 hs-rna.mat hs-rna.maf hs-rna.maf
    try last-genotype -v -S1M hs-rna.mat hs-rna.maf
    try last-genotype -b1 hs-rna.mat hs-rna.maf
    try last-genotype -b1 -c- hs-rna.mat hs-rna.maf
} 2>&1 |
grep -v version |
diff -u $(basename $0 .sh).out -
