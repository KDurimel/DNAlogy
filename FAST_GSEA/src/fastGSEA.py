#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, time, subprocess

__doc__="""
FastGSEA workfow. Performs GO-terms enrichment analysis between two gene sets. These gene sets
must be provided as files containing one international databank (ncbi, refseq, etc...) gene or
protein identifier per line.
@requires: U{python 2.7<https://www.python.org/downloads/>} or greater
@requires: U{Conda 4.4.10<https://conda.io/>} or greater with FastGSEA  conda environment (fastgsea.yml)
@requires: map.py
@requires: managefiles.sh
"""


def get_parser():
	"""
	Initlializing arguments.
	@return: arguments values
	@rtype: object
	"""

	parser = argparse.ArgumentParser(description='FastGSEA workfow. Performs GO-terms enrichment analysis between \
	two gene sets: sample and universe. These gene sets must be provided as files containing one international databank \
	(ncbi, refseq, etc...) gene or protein identifier per line. Supported ids are: \
	 \n \
	1. UniProtKB-AC \n \
	2. RefSeq \n \
	3. UniProtKB-ID \n \
	4. GeneID (EntrezGene) \n \
	5. GI \n \
	6. GO (as output only!) \n \
	7. UniRef100 \n \
	8. UniRef90 \n \
	9. UniRef50 \n \
	10. UniParc \n \
	11. UniGene \n \
	12. EMBL \n \
	13. EMBL-CDS \n \
	14. Ensembl \n \
	15. Ensembl_TRS \n \
	16. Ensembl_PRO')

	parser.add_argument('-ech', dest='ech', action="store",
						type=str, required=True, help='File containing the sample ids (one id per line) - REQUIRED')

	parser.add_argument('-univ', dest='univ', action="store",
						type=str, required=True, help='File containing the universe ids (one id per line) - REQUIRED')

	parser.add_argument('-mappingFile', dest='mappingFile', action="store",
						type=str, required=True, help='File used for offline mapping and multi ids support - REQUIRED')

	parser.add_argument('-fromOtherDB', dest='fromOtherDB', action="store",
						type=str, required=True, help='Ids supported by default are UniProtKB-AC and RefSeq ids. Use this \
						option if you want to use all the supported ids. Caution: results may be less reliable!')

	parser.add_argument('--mapOffline', dest='mapOffline', action='store_true', 
						help='Retrieve GO-terms without requesting Uniprot API\'s. \
						Faster but less reliable. Use this option in case of low internet bandwith, but \
						dont forget that the GO enrichment step (following the mapping step) will still \
						requires a reliable internet connection.', default=False)

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
        self._lifetime=time.time()

    @property
    def lifetime(self):
        return self._lifetime

    @lifetime.setter
    def lifetime(self, message_lifetime):
        self._lifetime  = message_lifetime + " " + str(int((time.time()-self._lifetime))) + " seconds"
        return self._lifetime 

# ************************************************************************************************************************************** #

def main():	
	################# Counter initialisation ########################
	compteur = Timer()

	##################### Get parser ###############################
	parser=get_parser()

	######## Print parser help if arguments missed #################
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	########### Manage workflow accorded to Args  ##################
	Arguments=parser.parse_args()

	########### Show time elapsed  ##################
	compteur.lifetime = "Script finished in"
	print compteur.lifetime


import map # custom libs: id mapping functions 

# map.any_ids_to_go(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
# map.ids_to_go(sys.argv[1],sys.argv[2],sys.argv[3])
# map.mk_susbet(sys.argv[1])
# map.ids_to_go_online(sys.argv[1],sys.argv[2],sys.argv[3])
# map.any_ids_to_go_online(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])



#MAIN

# if no subset file or if suset argument given : subset file

# if args = online : performs go retrieve online 

# args de 