#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function # python 3 like print
import os, sys, subprocess
import signal # used to handle gzip broken pipe when zgrep, more info : https://blog.nelhage.com/2010/02/a-very-subtle-bug/

__doc__="""
Id mapping script for a gene set. Takes as input a file (gene set) constitued by Refseq,Uniprot, or other supported ids.
Supported ids are:
1. UniProtKB-AC
2. UniProtKB-ID
3. GeneID (EntrezGene)
4. RefSeq
5. GI
6. PDB
7. GO
8. UniRef100
9. UniRef90
10. UniRef50
11. UniParc
12. PIR
13. NCBI-taxon
14. MIM
15. UniGene
16. PubMed
17. EMBL
18. EMBL-CDS
19. Ensembl
20. Ensembl_TRS
21. Ensembl_PRO
22. Additional PubMed

@requires: U{python 2.7<https://www.python.org/downloads/>} (tested with 2.7.6)
@requires: manageFiles.sh
"""

EXEC_PATH = sys.path[0] # Current running script path

def mk_susbet(idMappingFile):
	print("Map module: generating a subset from ids mapping file...",end="")
	cmd = 'zcat ' + idMappingFile + ' | cut -f 1,4,7 | gzip --stdout > ' + idMappingFile + '_subset.gz'
	subprocess.check_call(cmd, shell=True, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
	print("ok")


# Args : mode(offline ou online), file avec les ids
def ids_to_go():
	print("return to ids from refseq or uniprot ids.")

def ids_to_go_online():
	print("return to ids from refseq or uniprot ids.")

# whichdb : refseq, etc, outputprefix: go file output prefix
def any_ids_to_go(idMappingFile,idsFile,whichDb,outputPrefix):
	"""
	Make enrichment analysis for each comparison values
	@param compvalues: comparisons values returned by getcomp() 
	@type compvalues: string
	@param view: -view argument --> if view != 0 , graphical representation of enrichment will be processed
	@type view: int
	@param xmlname: the xml filename
	@type xmlname: string
	"""
	all_goterms = []
	with open(idsFile,'r') as ids:
		for id_ in ids:
			cmd = 'bash ' + EXEC_PATH+'/manageFiles.sh '+' --anythingToAnything '+ idMappingFile+' '+id_.strip()+' '+whichDb
			cmd2 = [EXEC_PATH+'/manageFiles.sh ',' --anythingToAnything ', idMappingFile+' ',id_.strip()+' ',whichDb]
			# print(cmd)
			#goterms = subprocess.check_output(cmd, shell=True, executable="/bin/bash", preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
			goterms = str(subprocess.check_output(cmd,shell=True,preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))).split(';')
			for goterm in goterms:
				goterm_clean = goterm.strip()
				# if goterms exist (not nonetype or something else bcs was processed by strip())
				if goterm_clean:
					all_goterms.append(goterm_clean)

	# write all goterms in a file
	with open (outputPrefix, "w+") as gofile:
		gofile.write('\n'.join(all_goterms))

def any_ids_to_go_online():
	print("blabla")
	# print "Map module: generating a subset from ids mapping file ...".strip()
	# bash_cmd = 'bash manageFiles.sh --anythingToAnything ../db/test.gz "Q6GZX1" "RefSeq"'
	# subprocess.check_output('bash manageFiles.sh --anythingToAnything ../db/test.gz "Q6GZX1" "RefSeq"',\
	# 	shell=True,preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL)).strip()
	# print "ok"