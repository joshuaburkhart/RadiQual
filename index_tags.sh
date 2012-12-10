#!/bin/bash

#Usage: ./rad_indexer.sh <unindexed rad file>

index=0;
while read line;
do
    echo "$index $line" >> rad_tags.idx.txt
    index=$((index+1))
    echo -n "."
done < $1
echo "done"
