#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail
#enter the workflow's final output directory ($1)
cd $1
#find all files, sort them, return their md5sums to std out
find . -name "*.vcf.gz" -xtype f | sort | while read -r f; do
    zcat "$f" | grep -v ^# | md5sum
done
ls | sed 's/.*\.//' | sort | uniq -c
