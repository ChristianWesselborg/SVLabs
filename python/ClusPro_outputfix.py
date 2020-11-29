# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 10:39:10 2020

@author: WESSEC

Format .pdb file outputs from cluspro to correctly load in models

purpose:
	line 1 reads "HEADER rec.pdb"
	line 1 rewritten as "HEADER rec.<model number>.pdb"
"""

import glob
import os

filedir = "c:/Users/WESSEC/Documents/Research/RPI GB/clusprooutput/6IRA_blocked/*/*/*.pdb"
dirlist = glob.glob(filedir)

for file in dirlist:
	print(file)
	
	f = open(file, 'r+')

	model = os.path.basename(f.name)
	modelnum = model.split(".")

	hold = f.readlines()
	if "rec.pdb" in hold[0]:
		hold[0] = hold[0].replace(".", "." + str(int(modelnum[-2])) + ".")
		f.seek(0)
		f.truncate()
		f.writelines(hold)
	f.close()

