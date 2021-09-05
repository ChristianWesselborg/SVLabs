# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 14:34:20 2020

@author: Christian

rough code for returning the residue list of a selection in PyMOL and saves to txt file in csv format
- used for getting list of all residues in a given selection/pbd file
"""
from pymol import cmd

def resadd(resn, resi, chain):
	try:
		if([resn, resi, chain] != reslist[-1]):
			reslist.append([resn,resi,chain])
	except:
		reslist.append([resn,resi,chain])

reslist = []
cmd.iterate("all","resadd(resn,resi,chain)") #first string is the selection you want to work with

hold = ""
for i in reslist:
	temp = str(i)
	temp = temp.translate({ord(x): None for x in "[]' "})
	temp = temp.split(",")
	hold += temp[0] + temp[1] + temp[2] + "\n"

f = open("c:/Users/WESSEC/Documents/Research/RPI GB/misc/pdb extract/3SQ9_reslist.txt", 'w')
f.write(hold)
f.close()
