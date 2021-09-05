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
cmd.load("c:/users/wessec/documents/research/rpi gb/pdbfiles/6ira.pdb")
cmd.alter("all", "b=0.0")

file = open("c:/users/wessec/documents/research/rpi gb/data/bfactor/6ira_blocked/zdock/mono_receptorWT.csv", 'r')
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
cmd.show_as("spheres", "interface")
cmd.color("grey", "not interface")

cmd.spectrum("b", "blue_red", "interface", minimum = 0, maximum = max)

cmd.select("top75", "b>" + str(max * 0.75))

cmd.label("top75 and n. CA", '"(%s%s%s)"%(resn,resi,chain)')
cmd.set("label_position",(3,0,20))

cmd.deselect()
