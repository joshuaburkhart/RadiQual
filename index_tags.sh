#!/bin/bash

#Usage: ./index_tags.sh <unindexed rad file>

index=0;
while read line;
do
    echo -e ">$index\n$line" >> rad_tags.idx.txt
    index=$((index+1))
    echo -n "."
done < $1
echo "done"
