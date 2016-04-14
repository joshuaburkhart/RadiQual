#!/usr/bin/env bash

#Usage: ./index_tags.sh <unindexed rad file>

if [ "$#" -ne 2 ]; then
    printf "ERROR: Too few arguments supplied. \
Please supply input and output filenames.\n\
Example usage: $ ./index_tags.sh ../data/input/rad_tags.txt ../data/input/rad_tags.idx.txt\n"
    exit
fi

index=0;
while read line;
do
    echo -e ">$index\n$line" >> $2
    index=$((index+1))
    echo -n "."
done < $1
echo "done"
