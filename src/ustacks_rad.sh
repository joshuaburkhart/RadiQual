#!/usr/bin/env bash

# For RAD files of the form:
# 0	0	0		0	+	consensus			TGCAGGGAAGCTATTAATTTAGACCTGGCGTCACGGTGTACAGGGCATG	0	0	1
# 0	0	0				model			OOOOOOOOOOOOOOOOOOEOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

if [ "$#" -ne 2 ]; then
    printf "ERROR: Too few arguments supplied. \
Please supply input and output filenames.\n\
Example usage: $ ./ustacks_rad.sh ../data/input/KC.tags.tsv ../data/input/rad_tags.txt\n"
    exit
fi

unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    cat $1 | \
    grep -P '^[0-9]+[\s]+[0-9]+[\s]+[0-9]+[\s]+[a-z]+[\s]+[0-9]+[\s]+[0-9a-zA-Z_-]+[\s]+[ACTG]+[\s]+' | \
    awk -F' ' '{print $NF}' | \
    uniq | \
    sort > $2
elif [[ "$unamestr" == 'Darwin' ]]; then
    cat $1 | \
    perl -nle 'print $& if m{^[0-9]+[\s]+[0-9]+[\s]+[0-9]+[\s]+[a-z]+[\s]+[0-9]+[\s]+[0-9a-zA-Z_-]+[\s]+[ACTG]+[\s]+}' | \
    awk -F' ' '{print $NF}' | \
    uniq | \
    sort > $2
else
    echo ERROR: "$unamestr" OS not supported
fi