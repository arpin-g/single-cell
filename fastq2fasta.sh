#!/bin/bash
cat $1 | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' | gzip -c - > $1.fasta.gz

