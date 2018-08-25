#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__="""
FastGSEA workfow. Performs GO-terms enrichment analysis between two gene sets. These gene sets
must be provided as files containing one international databank (ncbi, refseq, etc...) gene or
protein identifier per line.
@requires: U{python 2.7<https://www.python.org/downloads/>} or greater
@requires: U{Conda 4.4.10<https://conda.io/>} or greater with FastGSEA  conda environment (fastgsea.yml)
@requires: map.py
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
import map # homemade id mapping module : import all


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
						type=str, required=True, help='Sample ids file (one id per line): recquired')

	parser.add_argument('-univ', dest='univ', action="store",
						type=str, required=True, help='Universe ids file (one id per line): recquired')

	parser.add_argument('-mappingFile', dest='mappingFile', action="store",
						type=str, required=True, help='idmapping_selected.tab.gz file: recquired')

	parser.add_argument('-output', dest='output', action="store",
						type=str, required=True, help='output results prefix: recquired')

	parser.add_argument('-fromOtherDB', dest='fromOtherDB', action="store",
						type=str, required=False, help='Ids supported by default are UniProtKB-AC and RefSeq ids. Use this '+\
						'option if you want to use all the supported ids. Caution: results may be less reliable!')

	parser.add_argument('--mapOffline', dest='mapOffline', action='store_true', 
						help='Retrieve GO-terms without requesting Uniprot API\'s. '+\
						'Faster but less reliable. Use this option in case of low internet bandwith, but '+\
						'dont forget that the GO enrichment step (following the mapping step) will still '+\
						'requires a reliable internet connection.', default=False)

	parser.add_argument('--view', dest='view', action='store_true', 
						help='Generate graphical representation of the enrichment results (experimental feature)',
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
	SUPPORTED_IDS = ['UniProtKB-AC','RefSeq','UniProtKB-ID','GeneID',
	'GI','GO','UniRef100','UniRef90','UniRef50','UniParc','UniGene',
	'EMBL','EMBL-CDS','Ensembl','Ensembl_TRS','Ensembl_PRO']

	##################### Get parser ###############################
	parser=get_parser()

	######## Print parser help if arguments missed #################
	if len(sys.argv)==1:
		parser.print_help()
		exit(1)

	########### Manage workflow accorded to Args  ##################
	Arguments=parser.parse_args()

	################ Workflow starts here ##########################
	subprocess.check_call('mkdir -p ' + Arguments.output + '/tmp', shell=True)
	TMP_DIR = Arguments.output + '/tmp'
	if not exists(Arguments.mappingFile + '_subset.gz'):
		# Create idmapping subset file if it dont already exists
		map.mk_susbet(Arguments.mappingFile) 

	if Arguments.mapOffline:
		# Map ids files OFFLINE enabling only Refseq and GO ids support. Reliable solution
		if Arguments.fromOtherDB:
			if Arguments.fromOtherDB not in SUPPORTED_IDS:
				print 'fromOtherDB - bad argument: ' + Arguments.fromOtherDB + \
				'\nPlease use only supported ids. Program will stop now.'
				sys.exit(1)
			# Map ids files OFFLINE enabling all ids support. Results may be uncomplete
			else:
				map.any_ids_to_go(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt', Arguments.fromOtherDB) # Map sample ids
				map.any_ids_to_go(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt', Arguments.fromOtherDB) # Map universe ids
		# Map ids files OFFLINE enabling only Refseq and GO ids support. Reliable solution
		else:
			map.ids_to_go(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt') # Map sample ids
			map.ids_to_go(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt') # Map universe ids
	else:
		# Map ids files ONLINE enabling only Refseq and GO ids support. BEST solution for strong results
		if Arguments.fromOtherDB:
			if Arguments.fromOtherDB not in SUPPORTED_IDS:
				print 'fromOtherDB - bad argument: ' + str(Arguments.fromOtherDB) + \
				'\nPlease use only supported ids. Program will stop now.'
				sys.exit(1)
			else:
				map.any_ids_to_go_online(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt', Arguments.fromOtherDB) # Map sample ids
				map.any_ids_to_go_online(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt', Arguments.fromOtherDB) # Map universe ids
		# Map ids files ONLINE enabling all ids support. Results may be uncomplete
		else:
			map.ids_to_go_online(Arguments.mappingFile, Arguments.ech, TMP_DIR + '/go_ech_raw.txt') # Map sample ids
			map.ids_to_go_online(Arguments.mappingFile, Arguments.univ, TMP_DIR + '/go_univ_raw.txt') # Map universe ids

	# remove tmp files
	# subprocess.check_call('rm -r ' + Arguments.output + '/tmp', shell=True)
	################## Show time elapsed  ##########################
	TIMER.lifetime = "Workflow finished in"
	print TIMER.lifetime

# Executed only if script is not imported as a module
if __name__ == "__main__": 											
	main()											
