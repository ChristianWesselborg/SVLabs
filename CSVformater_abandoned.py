# -*- coding: utf-8 -*-
"""
The putsit in a CSV

Created on Wed Jul 22 12:30:55 2020

@author: Christian
"""

import os
import re

cdir = ""
masterlist = os.listdir()

sublist = list
for i in masterlist:
	if re.match("*output", masterlist[i]):
		sublist.append(masterlist[i])

for i in sublist:
	filer = open(os.path.normpath(cdir + "/" + sublist[i]), 'r')
	filew = open(os.path.normpath(cdir + "/datasheet.csv"), 'w')
	
	filer.readline()
	
	splitlist = filer.readline().strip('\n').split("    ")
	
	