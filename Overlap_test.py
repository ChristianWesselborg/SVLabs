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

#writepath = "c:/Users/WESSEC/Documents/Research/RPI GB/Data/interface overlap/zdock_6ira"
writepath = "c:/Users/WESSEC/Documents/Research/RPI GB/test2"

mainpath = "c:/Users/WESSEC/Documents/Research/RPI GB/zdockoutput/blocked_6ira/"

listsub = os.listdir(mainpath)
for folder in listsub:
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
				interface = interface[1].split(",")
				
				freqstore.append([int(modelnum)])
				
				selection = "model." + modelnum
				cmd.select(selection, "none")
				for model in interface[0:-1]:
					resn = model[0:3]
					resi = str(re.findall("\d+", model))[2:-2]
					chain = model[-1]
					
					cmd.select(selection, selection + " or (resn %s and resi %s and chain %s)" % (resn,resi,chain))
			line = data.readline()
		data.close()
		
		#select overlapping residues
		place = 0
		for modelnum in freqstore:
			cmd.select("overlap_model." + str(modelnum[0]), "selection and model." + str(modelnum[0]))
			
			overlap = cmd.get_model("overlap_model." + str(modelnum[0]))
			
			addfreq = 0
			for obj in range(len(overlap.atom)):
				if overlap.atom[obj-1].resi != overlap.atom[obj].resi: #is there a posibility of missing any residues here???
					addfreq += 1
			freqstore[place].append(addfreq)
			
			place += 1
		
		#write to file
		f = open(writepath + "/" + folder + ".txt", 'w')
		header = "original file: " + basefile
		f.write(header)
		header = "\nmodelnumber, #overlaping residues"
		f.write(header)
		for element in freqstore:
			f.write("\n" + str(element)[1:-1])
		f.close()
		
		#clean up the environment
		print("here")
		cmd.delete("all")
		cmd.fetch(receptormodel)
		select()

print("end")