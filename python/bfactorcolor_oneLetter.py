# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:52:00 2020

@author: Christian

rough code for coloring by b factor derived from csv file
"""

from pymol import cmd
import csv
import re

cmd.delete("all")

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

cmd.load("c:/users/wessec/documents/research/rpi gb/pdbfiles/2nao_sur/dim_stackedWT_sur.pdb")
cmd.alter("all", "b=0.0")

file = open("c:/users/wessec/documents/research/rpi gb/data/bfactor/dim_stackedligandWT.csv", 'r')
data = csv.reader(file)
next(data)

max = 0
for row in data:
	resi = int(re.findall(r'\d+', row[1])[0])
	resn = row[1][0:3]
	chain = row[1][-1]
	
	cmd.alter("resn %s and resi %s and chain %s"%(resn, resi, chain), "b=" + str(row[2]))
	if int(row[2]) > max:
		max = int(row[2])

print(max)

cmd.select("interface", selection = "b>0")
cmd.show_as("surface", "interface")
cmd.set("transparency", "0.3")
cmd.color("grey", "not interface")

cmd.spectrum("b", "blue_red", "interface", minimum = 0, maximum = max)

cmd.select("hetatms", "hetatm")
cmd.color("green", "hetatms")

cmd.select("top75", "b>" + str(max * 0.75))

cmd.label("top75 and n. CA", '"(%s)"%(one_letter[resn])')
cmd.set("label_position",(3,0,20))
