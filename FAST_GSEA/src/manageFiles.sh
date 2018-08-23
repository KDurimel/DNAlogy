#!/bin/bash -v
# Author: k.durimel@gmail.com
# Run bash commands lines for idmapping files management.

if [ -z "$1" ] || [ -z "$2" ] # $1 : processing workflow | $2 : file to process
  then
    echo "No sufficient args"
    exit 1
fi
case $1 in
  # Generates a file subset based only on Uniprot and refseq ids for faster data exploration
  "subset")
  zcat $2 | cut -f 1,4,7 | gzip --stdout > $2_subset.gz # Subset from ~6gb file
  ;;
  # Zgrep Uniprot or Refseq ids on the file subset (faster)
  "zgrep")
  if [ -z "$3" ]
    then
    echo "No sufficient args : pattern not found"
    exit 1
  fi
  zgrep -m 1 "$3" $2_subset.gz | cut -f 3 # Keep GO only if exists else null
  ;;
  # Zgrep whatever id to obtain GO-terms using the huge file
  "zgrep_whatever")
  zgrep -m 1 "$3" $2 | cut -f 7
esac
