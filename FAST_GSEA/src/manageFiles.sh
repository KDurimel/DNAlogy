#!/bin/bash -v
# Author: k.durimel@gmail.com
# Run bash commands lines for idmapping files management.

declare -A ids # declare ids as associative array (enables indexof and other specs)
ids=(['UniProtKB-AC']=1 ['UniProtKB-ID']=2 ['GeneID']=3 ['RefSeq']=4 ['GI']=5 ['PDB']=6 ['GO']=7 ['UniRef100']=8 ['UniRef90']=9 \
 ['UniRef50']=10 ['UniParc']=11 ['PIR']=12 ['NCBI-taxon']=13 ['MIM']=14 ['UniGene']=15 ['PubMed']=16 ['EMBL']=17 ['EMBL-CDS']=18 \
 ['Ensembl']=19 ['Ensembl_TRS']=20 ['Ensembl_PRO']=21 ['Additional Pubmed']=22)

# Please note users needs to generate a subset before perfoming --subset or --subsetToGo commands

if [ -z "$1" ] || [ -z "$2" ] # $1 : processing workflow | $2 : file to process
  then
    echo "No sufficient args"
    exit 1
fi
case $1 in
  # Generates a file subset based only on Uniprot and refseq ids for faster data exploration
  "--subset")
  zcat $2 | cut -f 1,4,7 | gzip --stdout > $2_subset.gz # Subset from ~6gb file
  ;;
  # Zgrep Uniprot or Refseq ids on the file subset (faster)
  "--subsetToGo")
  if [ -z "$3" ]
    then
    echo "No sufficient args : input id not given"
    exit 1
  fi
  zgrep -m 1 "$3" $2_subset.gz | cut -f 3 # Keep GO only if exists else null
  ;;
  # Zgrep whatever id to obtain GO-terms using the huge file
  "--anythingToGo")
  zgrep -m 1 "$3" $2 | cut -f 7
  # Zgrep whatever id to obtain GO-terms using the huge file
  ;;
  "--anythingToAnything")
  if [ -z "$4" ]
    then
    "No sufficient args : output id not given"
    exit 1
  fi
  zgrep -m 1 "$3" $2 | cut -f ${ids[$4]} #$3 input ids to grep ; $2 mode ; $4 output ids
  ;;
esac
