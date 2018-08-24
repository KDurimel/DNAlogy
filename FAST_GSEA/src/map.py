#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function # python 3 like print
import os, sys, time, subprocess
import signal # Async alerts module used to handle gzip broken pipe when zgrep, more info : https://blog.nelhage.com/2010/02/a-very-subtle-bug/
import urllib,urllib2,requests	# HTML requests
import re #regex

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
	print('Map module: generating a subset from ids mapping file...', end = '')
	cmd = 'zcat ' + idMappingFile + ' | cut -f 1,4,7 | gzip --stdout > ' + idMappingFile + '_subset.gz'
	subprocess.check_call(cmd, shell=True, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))
	print('ok')


def ids_to_go(idMappingFile, idsFile, outputPrefix):
	"""
	Takes as input Uniprot or Refseq ids and returns its corresponding GO-terms. Please use this 
	function if you want to retrieve GO-terms using a input file containing Refseq or Uniprot ids.
	The id mapping file subset is automatically generated if it was not already the case (e.g First
	script launch , or file deleted.) 
	@param idMappingFile: path to the .gz file used for id mapping 
	@type idMappingFile: string
	@param idsFile: path to the input file containing all the ids we have to map 
	@type view: string
	@param outputPrefix: file path output prefix
	@type outputPrefix: string
	"""
	print('Process: offline GO mapping from RefSeq or Uniprot IDs...', end = '')
	all_goterms = [] # This will contains all the go-terms retrieved
	with open(idsFile,'r') as ids:
		for id_ in ids:
			
			cmd = 'bash ' + EXEC_PATH + '/manageFiles.sh '+ ' --subsetToGo ' + idMappingFile + ' ' + id_.strip()
			goterms = str(subprocess.check_output(cmd, shell=True, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))).split(';')
			for goterm in goterms:
				
				goterm_clean = goterm.strip()
				# if goterm exist (not nonetype or something else bcs was processed by strip())
				if goterm_clean:
					all_goterms.append(goterm_clean)

	# write all goterms in a file
	with open (outputPrefix, 'w') as gofile:
		gofile.write('\n'.join(all_goterms))

	print('ok')

def ids_to_go_online(idMappingFile, idsFile, outputPrefix):
	"""
	This function follows same principles as ids_to_go() but performs id mapping online. It the best
	choice for reliability and up-to-date GO-terms. Ids_to_to_online() always requires idMappingFile
	parameter to retrieve refseq ids offline. 
	Takes as input Uniprot or Refseq ids and returns its corresponding GO-terms. Please use this 
	function if you want to retrieve GO-terms using a input file containing Refseq or Uniprot ids.
	@param idMappingFile: path to the .gz file used for id mapping 
	@type idMappingFile: string
	@param idsFile: path to the input file containing all the ids we have to map 
	@type view: string
	@param outputPrefix: file path output prefix
	@type outputPrefix: string
	"""
	print('Process: Online GO mapping from RefSeq or Uniprot IDs...', end = '')
	all_goterms = [] # This will contains all the go-terms retrieved

	with open(idsFile,'r') as ids:
		for id_ in ids:
			
			#if this line contains a Refseq id (NP_xxx, Wp_xxx...) performs offline id mapping
			if 'P_' in id_.strip():
				print(id_)
				cmd = 'bash ' + EXEC_PATH + '/manageFiles.sh ' + ' --anythingToGo ' + idMappingFile + ' ' + id_.strip()
				goterms = str(subprocess.check_output(cmd,shell=True,preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))).split(';')
				for goterm in goterms:
					
					goterm_clean = goterm.strip()
					# if goterm exist (not nonetype or something else bcs was processed by strip())
					if goterm_clean:
						all_goterms.append(goterm_clean)

			# If this line contains a Uniprot ID, perform online id mapping
			else:
				dicogo = {} # used to avoid GO-terms duplicates for a same Uniprot id (it really happens sometimes)
				res_text = None # True if GET request received a response
				while res_text == None:

					try:
						url='https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=' + id_
						res = requests.get(url, headers = { "Accept" : "text/tsv"}, verify = True)
						res_text = res.text
						goterms = re.findall(r'\bGO:\d+',res_text) # re.findall returns a list of GO ids strings
						
						for goterm in goterms:

							if goterm not in dicogo.keys() and goterm!='GO':
								all_goterms.append(goterm)
								dicogo[goterm]=id_
								print(goterm)

					except urllib2.URLError:
						print('Timeout, we\'ll wait 5 seconds before re-trying the request...')
						time.sleep(5)
						pass	

				dicogo[goterm]= [] # cleared for next loop
	# write all goterms in a file
	with open (outputPrefix, 'w') as gofile:
		gofile.write('\n'.join(all_goterms))
	print('ok')



def any_ids_to_go(idMappingFile, idsFile, outputPrefix, whichDb):
	"""
	Takes as input any suppported ids and returns its corresponding GO-terms. Slower
	than ids_go_to beacause it do not use the idMappingFile subset and requests id per id.
	@param idMappingFile: path to the .gz file used for id mapping 
	@type idMappingFile: string
	@param idsFile: path to the input file containing all the ids we have to map 
	@type view: string
	@param whichDb: the ids database prefix that we want in output (for example, if we 
	want go ids, whichDB = GO. Please refer to "supported ids" at the top section)
	@type whichDb: string
	@param outputPrefix: file path output prefix
	@type outputPrefix: string
	"""
	print('Process: Offline GO mapping from any IDs...', end = "")
	all_goterms = [] # This will contains all the go-terms retrieved

	with open(idsFile,'r') as ids:
		for id_ in ids:

			cmd = 'bash ' + EXEC_PATH + '/manageFiles.sh '+ ' --anythingToAnything '+ idMappingFile + ' ' + id_.strip() + ' ' + whichDb
			goterms = str(subprocess.check_output(cmd, shell = True, preexec_fn = lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))).split(';')
			for goterm in goterms:

				goterm_clean = goterm.strip()
				# if goterms exist (not nonetype or something else bcs was processed by strip())
				if goterm_clean:
					all_goterms.append(goterm_clean)

	# write all goterms in a file
	with open (outputPrefix, 'w') as gofile:
		gofile.write('\n'.join(all_goterms))
	print('ok')

def any_ids_to_go_online(idMappingFile, idsFile, outputPrefix, whichDb):
	"""
	This function follows same principles as any_ids_to_go() but performs id mapping online.
	This function is an experimental feature and support only the ids supported by the 
	U{ Uniprot API <https://www.ebi.ac.uk/QuickGO/api/index.html>} 
	Takes as input any suppported ids and returns its corresponding GO-terms. Slower
	than any_ids_go_to() beacause it do not use the idMappingFile subset and requests id per id.
	@param idMappingFile: path to the .gz file used for id mapping 
	@type idMappingFile: string
	@param idsFile: path to the input file containing all the ids we have to map 
	@type view: string
	@param whichDb: the ids database prefix that we want in output (for example, if we 
	want go ids, whichDB = GO. Please refer to "supported ids" at the top section)
	@type whichDb: string
	@param outputPrefix: file path output prefix
	@type outputPrefix: string
	"""
	print('Process: Online GO mapping from any IDs (experimental feature, use at your own risks!) ...', end = "")
	all_goterms = [] # This will contains all the go-terms retrieved

	with open(idsFile,'r') as ids:
		for id_ in ids:
			
			#if this line contains a Refseq id (NP_xxx, Wp_xxx...) performs offline id mapping
			if 'P_' in id_.strip():
				print(id_)
				cmd = 'bash ' + EXEC_PATH + '/manageFiles.sh ' + ' --anythingToAnything ' + idMappingFile + ' ' + id_.strip() + ' ' + whichDb
				goterms = str(subprocess.check_output(cmd,shell=True,preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))).split(';')
				for goterm in goterms:
					
					goterm_clean = goterm.strip()
					# if goterm exist (not nonetype or something else bcs was processed by strip())
					if goterm_clean:
						all_goterms.append(goterm_clean)

			# If this line contains a Uniprot ID, perform online id mapping
			else:
				dicogo = {} # used to avoid GO-terms duplicates for a same Uniprot id (it really happens sometimes)
				res_text = None # True if GET request received a response
				while res_text == None:

					try:
						url='https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=' + id_
						res = requests.get(url, headers = { "Accept" : "text/tsv"}, verify = True)
						res_text = res.text
						goterms = re.findall(r'\bGO:\d+',res_text) # re.findall returns a list of GO ids strings
						
						for goterm in goterms:

							if goterm not in dicogo.keys() and goterm!='GO':
								all_goterms.append(goterm)
								dicogo[goterm]=id_
								print(goterm)

					except urllib2.URLError:
						print('Timeout, we\'ll wait 5 seconds before re-trying the request...')
						time.sleep(5)
						pass	

				dicogo[goterm]= [] # cleared for next loop
	# write all goterms in a file
	with open (outputPrefix, 'w') as gofile:
		gofile.write('\n'.join(all_goterms))
	print('ok')