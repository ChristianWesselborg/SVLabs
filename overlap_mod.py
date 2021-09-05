# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 13:17:46 2021

@author: WESSEC
"""


# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:28:25 2021

@author: WESSEC

Purpose: Load proximal residue lists from active compounds search to appropriate model
"""

#ititialize
from pymol import cmd

import os
#import glob
#import re

#clear PyMOL
cmd.delete("all")

#set path where comparison selection is found
compath = "c:/Users/wessec/documents/research/rpi gb/misc/proximal residue lists/4pe5/"

filelist = os.listdir(compath) #stores list of files in compath

#load correct model

model = os.path.basename(os.path.normpath(compath))
print(model)

cmd.fetch(model)

#run through filelist saving in selectstore

selectstore = []

switcher = {
	"A":"C",
	"B":"D",
	"C":"A",
	"D":"B"
	}

for file in filelist:
	current = open(compath + file)
	
	line = current.readline()
	while line:
		if line[:-1] in selectstore:
			pass
		else:
			chain = line[-2]
			
			identicalchain = switcher.get(chain, "Error")
			selectstore.append(line[:-1])
			
			fliped = line[:-2] + identicalchain
			selectstore.append(fliped)
		line = current.readline()
	current.close()


#select every entry in selectstore
cmd.select("selection", "none")

for element in selectstore:
	name = element[0:3]
	id = element[3:-1]
	
	cmd.select("selection", "selection or resn %s and resi %s" % (name, id))