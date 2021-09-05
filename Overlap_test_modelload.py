# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:28:25 2021

@author: WESSEC

Purpose: Identify overlap of one selection and save results as a text file ranking each model by overlap
"""

#ititialize
from pymol import cmd

import os
import glob
import re

#clear PyMOL
cmd.delete("all")

#set path where comparison selection is found
compath = "c:/Users/wessec/documents/research/rpi gb/misc/proximal residue lists/6ira/"

filelist = os.listdir(compath) #stores list of files in compath

#load correct model

receptormodel = os.path.basename(os.path.normpath(compath))
print(receptormodel)

cmd.fetch(receptormodel)

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
def select ():
	cmd.select("selection", "none")
	
	for element in selectstore:
		name = element[0:3]
		id = element[3:-1]
		
		cmd.select("selection", "selection or resn %s and resi %s" % (name, id))


select()



##find overlap between model and selection
#store all data directories

mainpath = "c:/Users/WESSEC/Documents/Research/RPI GB/zdockoutput/blocked_6ira/"

onelist = ["dim_stackedA2T"]
for folder in onelist:
	outputlist = glob.glob(mainpath + folder + "/*_output.txt")
	print(folder)
	
	#iterate through all data directories
	for output in outputlist:
		basefile = os.path.basename(output)
		data = open(output)
		
		#store all valid models in an array to have overlaps attached
		freqstore = []
		
		line = data.readline() #skips formating information in file
		line = data.readline() #starts at first line of data
		while line:
			if "invalid" not in line:
				interface = line.split("    ")
				
				modelnum = interface[0]
				
				cmd.load(mainpath + folder + "/predictions/" + folder + "." + modelnum + ".pdb")
			line = data.readline()
		data.close()
		
		#select overlapping residues

		
		#write to file


print("end")