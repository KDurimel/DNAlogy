#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, subprocess

__doc__="""
FastGSEA workfow. Performs GO-terms enrichment analysis between two gene sets.

@requires: U{python 2.7<https://www.python.org/downloads/>} (tested with 2.7.6)
@requires: map.py
@requires: managefiles.sh
@requires: xxx
@requires: xxx
"""

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