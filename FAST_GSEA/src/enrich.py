#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__="""
GSEA script. Launches Enrich.R script and manage its results.
@requires: U{python 2.7<https://www.python.org/downloads/>} or greater
@requires: U{R 3.3.3 RC <https://www.r-project.org/>}
@requires: U{Conda 4.4.10<https://conda.io/>} or greater with FastGSEA  conda environment (fastgsea.yml)
"""

__author__ = 'Kevin Durimel'
__credits__ = ['Kevin La', 'Arnaud Felten', 'Meryl Vila-Nova', 'Nicolas Radomski']
__license__ = 'GPL'
__version__ = '3.0'
__maintainer__ = 'Kevin Durimel'
__email__ = 'k@durimel@gmail.com'
__status__ = 'Alpha'

import sys, subprocess # too many uses --> import all.

EXEC_DIR = sys.path[0] # Current running script path

def launchGSEA(outputDir):
	print 'Online GO enrichment and hypergeometric tests:'
	subprocess.check_call('R --vanilla --slave --args ' + outputDir + ' < ' + EXEC_DIR + "/enrich.R", shell = True)
	# Add header
	subprocess.check_call(' echo "GO:ID;Go term;Number of hits;Expected number of hits;Go level;P-value;Corrected p-value;Aspect" > ' + \
		outputDir + '/../hyperesults.csv', shell = True)
	# Merge all files in one final file
	subprocess.check_call('cat ' + outputDir + '/hyperesults_* >> ' +  outputDir + '/../hyperesults.csv', shell = True) 
	print '...ok'