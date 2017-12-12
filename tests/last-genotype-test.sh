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
    try last-genotype -M0 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -m5 hs-rna.mat hs-rna.maf
    try last-genotype -p1 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -p3 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -f1e6 -s50 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -f1000 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -p'3,chr1*:1' hs-rna.mat hs-rna.maf
    try last-genotype -M0 -f1000 hs-rna.mat hs-rna.maf hs-rna.maf
    try last-genotype -M0 -v -S1M hs-rna.mat hs-rna.maf
    try last-genotype -M0 -b1 hs-rna.mat hs-rna.maf
    try last-genotype -M0 -b1 -c- hs-rna.mat hs-rna.maf
} 2>&1 |
grep -v version |
diff -u $(basename $0 .sh).out -
