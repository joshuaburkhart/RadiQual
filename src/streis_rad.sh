#!/usr/bin/env bash

# For RAD files of the form:
# Catalog ID	Consensus Sequence	Num Parents	Num SNPs	SNPs	Num Alleles	Alleles	Deleveraged	ucsd1	ucsd2	ucsd6	ucsd7	ucsd8	elf1	elf	elf6	elf7	elf8	elf9	lh1	lh2	lh3	lh4	lh5	lh6	lo1	lo2	lo3	lo4	lo5	lo6	pct1	pct2	pct3	pct4	pct8	pct	potr1	potr2	potr3	potr4	potr8	potr9
# 1	TGCAGAGGAAAAGCTTTATCCCAGCGAATAGAGGATGCTCAGCTCCGCACAACATAGGGATTGGCACACTTTCGGGAAAAACATAGAATAGTATC	35	0		1		1	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus	consensus

if [ "$#" -ne 2 ]; then
    printf "ERROR: Too few arguments supplied. \
Please supply input and output filenames.\n\
Example usage: $ ./streis_rad.sh ../data/input/raw_stacks_6.tsv.no_crs ../data/input/rad_tags.txt\n"
    exit
fi

unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    cat $1 | \
    grep -P '^[0-9]+[\s]+[ACTG]+[\s]+[0-9]+[\s]+[0-9]+[\s]+1[\s]+' | \
    awk -F' ' '{print ">"$1"\n"$2}' > $2
elif [[ "$unamestr" == 'Darwin' ]]; then
    cat $1 | \
    perl -nle 'print $& if m{^[0-9]+[\s]+[ACTG]+[\s]+[0-9]+[\s]+[0-9]+[\s]+1[\s]+}' | \
    awk -F' ' '{print ">"$1"\n"$2}' > $2
else
    echo ERROR: "$unamestr" OS not supported
fi