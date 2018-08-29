#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__="""
FastGSEA workfow. Performs GO-terms enrichment analysis between two gene sets. These gene sets
must be provided as files containing one international databank (ncbi, refseq, etc...) gene or
protein identifier per line.
@requires: U{python 2.7<https://www.python.org/downloads/>} or greater
@requires: U{R 3.3.3 RC <https://www.r-project.org/>}
@requires: U{Conda 4.4.10<https://conda.io/>} or greater with FastGSEA  conda environment (fastgsea.yml)
@requires: map.py
@requires: enrich.py
@requires: trim.py
@requires: managefiles.sh
"""

__author__ = 'Kevin Durimel'
__credits__ = ['Kevin La', 'Arnaud Felten', 'Meryl Vila-Nova', 'Nicolas Radomski']
__license__ = 'GPL'
__version__ = '3.0'
__maintainer__ = 'Kevin Durimel'
__email__ = 'k@durimel@gmail.com'
__status__ = 'Alpha'

######## IMPORT BUILT IN MODULES ########
from os.path import exists # verify if path exists
from time import time # handle system time
from argparse import RawDescriptionHelpFormatter, ArgumentParser # help text formatter, argument parser
import sys, subprocess # too many uses --> import all.

######## IMPORT HOMEMADE MODULES ########
import map # homemade id mapping module : import all functionalities
from enrich import launchGSEA # homemade go-enrichment module
import trim # homemade ontology files management : import all functionalities


def get_parser():
	"""
	Initlializing arguments.
	@return: arguments values
	@rtype: object
	"""

	parser = ArgumentParser(description='FastGSEA workfow. Performs GO-terms enrichment analysis between\
	two gene sets: sample and universe. These gene sets must be provided as files containing one international databank\
	(ncbi, refseq, etc...) gene or protein identifier per line. Supported ids are:\n\
	1. UniProtKB-AC \n\
	2. RefSeq \n\
	3. UniProtKB-ID \n\
	4. GeneID (EntrezGene) \n\
	5. GI \n\
	6. GO (as output only!) \n\
	7. UniRef100 \n\
	8. UniRef90 \n\
	9. UniRef50 \n\
	10. UniParc \n\
	11. UniGene \n\
	12. EMBL \n\
	13. EMBL-CDS \n\
	14. Ensembl \n\
	15. Ensembl_TRS \n\
	16. Ensembl_PRO',
	formatter_class=RawDescriptionHelpFormatter)

	parser.add_argument('-ech', dest='ech', action="store",
						type=str, required=True, help='Sample ids file (one id per line): required')

	parser.add_argument('-univ', dest='univ', action="store",
						type=str, required=True, help='Universe ids file (one id per line): required')

	parser.add_argument('-mappingFile', dest='mappingFile', action="store",
						type=str, required=True, help='idmapping_selected.tab.gz file: required')

	parser.add_argument('-output', dest='output', action="store",
						type=str, required=True, help='output results prefix: required')

	parser.add_argument('--fromOtherDB', dest='fromOtherDB', action='store_true', 
						help='Ids supported by default are UniProtKB-AC and RefSeq ids. Use this '+\
						'option if you want to activate all ids support (slower). Caution: results may be less reliable!',
						default=False)

	parser.add_argument('-obo', dest='obo', action="store",
						type=str, required=False, help='Gene ontology .obo graph file used when "--trim" option activated'+\
						'(old go-basic.obo or latest gosubset_prok.obo)')

	parser.add_argument('-toDB', dest='toDB', action="store",
						type=str, required=False, help='databank identifier wanted (e.g GO,Uniref100...)'  +\
						'as output when "--mapOnly" option activated')

	parser.add_argument('-fastmode', dest='fastmode', action="store",
						type=str, required=False, help='optimize fastGOEA speeds for the entered identifier' +\
						'(Refseq,UniProtKB-AC,etc...)\nthis option improves execution speeds but disables multi-ids support')

	parser.add_argument('--mapOffline', dest='mapOffline', action='store_true', 
						help='Retrieve GO-terms without requesting Uniprot API\'s. '+\
						'Faster but less reliable. Use this option in case of low internet bandwith, but '+\
						'dont forget that the GO enrichment step (following the mapping step) will still '+\
						'requires a reliable internet connection.', default=False)

	parser.add_argument('--trim', dest='trim', action='store_true', 
						help='Trim prokaryotic GO-terms. This feature may be DEPRECATED in the future! More information here: '+\
						'https://github.com/geneontology/go-ontology/issues/16077',
						default=False)

	parser.add_argument('--view', dest='view', action='store_true', 
						help='Outputs the enrichment results in a 2D graph (experimental feature)',
						default=False)

	parser.add_argument('--mapOnly', dest='mapOnly', action='store_true', 
						help='Perform only the "MAP" step and keep its results',
						default=False)

	parser.add_argument('--keepTmp', dest='keepTmp', action='store_true', 
						help='Keep temporary files folder(tmp).',
						default=False)


	return parser


class Timer(object):
	"""
	Initlialize and manage timers.
	@return: started time or time elapsed
	@rtype: object
	"""
	def __init__(self):
		self._lifetime=time()

	@property
	def lifetime(self):
		return self._lifetime

	@lifetime.setter
	def lifetime(self, message_lifetime):
		self._lifetime  = message_lifetime + " " + str(int((time()-self._lifetime))) + " seconds"
		return self._lifetime


# *************************************************************************************************************************************

def main():	
	############# Global varuables initialisation ##################
	TIMER = Timer()
	EXEC_DIR = sys.path[0] # Current running script path
	# idmapping gz columns (used for retrieve cut column when mksubset)
	ALL_IDS = ['UniProtKB-AC','UniProtKB-ID','GeneID','RefSeq',
	'GI','PDB','GO','UniRef100','UniRef90','UniRef50','UniParc',
	'PIR','NCBI-taxon','MIM','UniGene','Pubmed','EMBL','EMBL-CDS','Ensembl',
	'Ensembl_TRS','Ensembl_PRO','Additional Pubmed']
	# supported ids : used to control user args input
	SUPPORTED_IDS = ['UniProtKB-AC','RefSeq','UniProtKB-ID','GeneID',
	'GI','GO','UniRef100','UniRef90','UniRef50','UniParc','UniGene',
	'EMBL','EMBL-CDS','Ensembl','Ensembl_TRS','Ensembl_PRO']

	##################### Get parser ###############################
	parser=get_parser()

	######## Print parser help if arguments missed #################
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	########### Manage workflow accorded to Args  ##################
	Arguments=parser.parse_args()

	################ Workflow starts here ##########################

	# 	************************ 
	# 	****** ID MAPPING ****** 
	# 	************************ 
	subprocess.check_call('mkdir -p ' + Arguments.output + '/tmp', shell=True)
	TMP_DIR = Arguments.output + '/tmp'

	#
	# IF FASTMODE ACTIVATED : Generate subset and use it for mapping (consequence :  no multi-ids support)
	if Arguments.fastmode:
		# Create idmapping subset if supported and its not GO (bcs why a subset for mapping go to go?)
		if Arguments.fastmode in SUPPORTED_IDS:
			# if used want to create a subset for GO ids its a mistake
			if Arguments.fastmode == "GO":
				parser.print_help()
				print '\nERROR :  generate subset for go ids is not allowed' + \
				'\nPlease use only supported ids. Program will stop now.'
				sys.exit(1)
			else:
				map.mk_susbet(Arguments.mappingFile,int(ALL_IDS.index(Arguments.fastmode))+1,Arguments.fastmode) 

	#
	# MAPPING ONLY MODE
	#
	if Arguments.mapOnly:
		if Arguments.toDB:
			if Arguments.toDB in SUPPORTED_IDS:
				map.any_ids_to_any_ids(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/../ech_mapped_ids.txt', Arguments.toDB) # Map sample ids
				map.any_ids_to_any_ids(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/../univ_mapped_ids.txt', Arguments.toDB) # Map universe ids
				# Remove tmp files
				if not Arguments.keepTmp:
					subprocess.check_call('rm -r ' + TMP_DIR, shell = True)
				# Map only : so exit without error
				sys.exit(0)
			else:
				parser.print_help()
				print '\nERROR : toDB - bad argument: ' + Arguments.toDB + \
				'\nPlease use only supported ids. Program will stop now.'
				sys.exit(1)
		else:
			parser.print_help()
			print '\nERROR : -toDB argument missed'
			sys.exit(1)

	#
	# FULL MODE : ONLINE OR OFFLINE
	#

	ECH_OUTPUT = 'go_ech_raw.txt' # ECH_OUTPUT = trimmed or not trimmed reads (values changes during the workflow if --trim option was selected)
	UNIV_OUTPUT = 'go_univ_raw.txt'  # UNIV_OUTPUT = trimmed or not trimmed reads (values changes during the workflow if --trim option was selected)

	if Arguments.mapOffline:
		# Map ids files OFFLINE enabling only all ids support.
		if Arguments.fromOtherDB:
			map.any_ids_to_go(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/' + ECH_OUTPUT) # Map sample ids
			map.any_ids_to_go(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/' + UNIV_OUTPUT) # Map universe ids
		# Map ids files OFFLINE enabling only Refseq and GO ids support. Faster and most reliable solution
		else:
			map.ids_to_go(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt') # Map sample ids
			map.ids_to_go(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt') # Map universe ids
	else:
		# Map ids files ONLINE enabling all ids support. Results may be uncomplete
		if Arguments.fromOtherDB:
			map.any_ids_to_go_online(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt') # Map sample ids
			map.any_ids_to_go_online(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt') # Map universe ids
		# Map ids files ONLINE enabling only Refset and Uniprot ids support. BEST solution for strong results
		else:
			map.ids_to_go_online(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt') # Map sample ids
			map.ids_to_go_online(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt') # Map universe ids


	# Trim prokarytic GO-terms if asked by user
	if Arguments.trim:
		if Arguments.obo:
			# Automatically watch if a subset file exists, and generates it if its not the case
			trim.mk_subset(Arguments.obo, TMP_DIR + '/gosubset.txt') 
			# Trim non prokaryote and non obsolete terms from enrichment results
			trim.trim(TMP_DIR + '/gosubset.txt', TMP_DIR + '/go_ech_raw.txt')
			trim.trim(TMP_DIR + '/gosubset.txt', TMP_DIR + '/go_univ_raw.txt')
			# Set new file input for launchGSEA (enrichment)
			ECH_OUTPUT = 'go_ech_raw.txt_cleaned.txt'
			UNIV_OUTPUT = 'go_univ_raw.txt_cleaned.txt'
		else:
			parser.print_help()
			print "\nERROR : Please provide a obo file!"
			exit(1)

	# 	************************ 
	# 	**** GO ENRICHMENT ***** 
	# 	************************ 
	# Gene set enrichment and hypergeometric tests using R scripts called by python map module
	launchGSEA(TMP_DIR,ECH_OUTPUT,UNIV_OUTPUT)
	
	# Generates GO distribution plot
	if Arguments.view:
		subprocess.check_call('R --vanilla --slave --args ' + TMP_DIR + '/../hyperesults.csv < ' + EXEC_DIR + '/goView.R', shell = True)
	# Remove tmp files
	if not Arguments.keepTmp:
		subprocess.check_call('rm -r ' + TMP_DIR, shell = True)

	################## Show time elapsed  ##########################
	TIMER.lifetime = "Workflow finished in"
	print TIMER.lifetime

# Executed only if script is not imported as a module
if __name__ == "__main__": 											
	main()											
