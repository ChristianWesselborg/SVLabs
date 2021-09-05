'''
Program loads a list of pdb files from the current folder with name suffix.*.pdb, where suffix is the common 'suffix' for all complexes, and * is an integer specified by min and max.

USAGE: Open a pymol window, and in the commandline type: run load_complexes.py

INPUT PARAMETERS: (can be assigned below in section parameter initialization)
 (1) suffix (common for all complexes) 
 (2) chrec (the chains for the receptor) 
 (3) chlig (the chains for the ligand) 
 (4) min, max (lower and upper limits for predictions to be loaded) 
 (5) intcut (Distance cutoff for what should be considered interface (in Angstrom))

OUTPUT: A text file - intresidues.txt
 FORMAT: 
  Each line of the file contains the complex number, the receptor interface residues list, the ligand interface residues list. 
  The residue list is comma separated and each residue is designated by (ResidueName)(ResidueNumber)(ResidueChain)

RESTRICTIONS/LIMITATION: 
Assumes all necessary input files (docking results)  are present in the current directory, including this script
Assumes the chain ids of receptor and ligand are unique (no common chain IDs between receptor and ligand) 
'''

from pymol import cmd

import re

#GUI setup
import os
import tkinter
import tkinter.font as tkFont
from tkinter import filedialog
from datetime import datetime

cwd = os.getcwd()
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")

#save settings function
def write():
	main.destroy()
	if(save.get() == 0):
		try:
			if(os.path.exists(os.path.normpath(g.get())) == True and g.get() != ""):
				filew = open(g.get() + "/" + now + ".txt", 'w')
			else:
				if(os.path.exists(os.path.normpath(h.get())) == True and h.get() != ""):
					filew = open(h.get() + "/" + now + ".txt", 'w')
				else:
					raise OSError
			filew.write("time and date: (YYYYMMDD_HHMMSS) " + now + "\n")
			filew.write("file prefix>" + a.get() + "\n")
			filew.write("receptor chains>" + b.get() + "\n")
			filew.write("ligand chains>" + c.get() + "\n")
			filew.write("file range>" + d.get() + "\n")
			filew.write("excluded files>" + e.get() + "\n")
			filew.write("interaction radius>" + f.get() + "\n")
			filew.write("results directory>" + g.get() + "\n")
			filew.write(".pbd file directory>" + h.get() + "\n")
			filew.close()
		except Exception:
			er = tkinter.Tk()
			er.geometry("400x300+%d+%d" % (x, y))
			message = tkinter.Label(er, text = "Error in saving settings file. Please check that the correct directory was entered", wraplength = 350, justify = "center")
			message.pack()
			traceres = tkinter.Label(er, text = "results directory tried: " + g.get())
			traceres.pack()
			tracewkdir = tkinter.Label(er, text = ".pdb directory tried: " + h.get())
			tracewkdir.pack()
			close = tkinter.Button(er, text = "Close", command = er.destroy)
			close.pack()
	complete.set(1)

def iset():
	main.choosefile = filedialog.askopenfilename(title = "Select settings file to load", filetypes = [("Text files", "*.txt")])
	i.set(main.choosefile)

def choosefile(passed):
	main.choosefile = filedialog.askdirectory(title = "Select directory")
	if(passed == 0):
		g.set(main.choosefile)
	else:
		h.set(main.choosefile)

def load():
	try:
		filer = open(i.get(), 'r')
		filer.readline()
		a.set(filer.readline().split(">")[1].strip('\n'))
		b.set(filer.readline().split(">")[1].strip('\n'))
		c.set(filer.readline().split(">")[1].strip('\n'))
		d.set(filer.readline().split(">")[1].strip('\n'))
		e.set(filer.readline().split(">")[1].strip('\n'))
		f.set(filer.readline().split(">")[1].strip('\n'))
		g.set(filer.readline().split(">")[1].strip('\n'))
		h.set(filer.readline().split(">")[1].strip('\n'))
	except Exception:
		er = tkinter.Tk()
		er.geometry("400x300+%d+%d" % (x, y))
		message = tkinter.Label(er, text = "Error in reading settings file. Please ensure integrity of loaded settings file", wraplength = 350, justify = "center")
		message.pack()
		close = tkinter.Button(er, text = "Close", command = er.destroy)
		close.pack()

def clr():
	a.set("")
	b.set("")
	c.set("")
	d.set("")
	e.set("")
	f.set("")
	g.set("")
	h.set("")
	i.set("")

#GUI interface
main = tkinter.Tk()
main.title("Settings configuration")

w = 700
h = 400
ws = main.winfo_screenwidth() # width of the screen
hs = main.winfo_screenheight() # height of the screen
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)
main.geometry("%dx%d+%d+%d" % (w, h, x, y))

fontset = tkFont.Font(size = 12)

a = tkinter.StringVar(main)
b = tkinter.StringVar(main)
c = tkinter.StringVar(main)
d = tkinter.StringVar(main)
e = tkinter.StringVar(main)
f = tkinter.StringVar(main)
g = tkinter.StringVar(main)
h = tkinter.StringVar(main)
i = tkinter.StringVar(main)
save = tkinter.IntVar(main)
complete = tkinter.IntVar(main, value = 0)

sufprompt = tkinter.Label(main, text = "common file prefix")
sufprompt.configure(font = fontset)
sufprompt.grid(row = 0, column = 0, columnspan = 2, sticky = "W")
sufent = tkinter.Entry(main, bd = 5, textvariable = a)
sufent.configure(font = fontset)
sufent.grid(row = 0, column = 1, sticky = "WE")

recprompt = tkinter.Label(main, text = "receptor chains ex: \"A+B+C...\"")
recprompt.configure(font = fontset)
recprompt.grid(row = 1, column = 0, columnspan = 2, sticky = "W")
recent = tkinter.Entry(main, bd = 5, textvariable = b)
recent.configure(font = fontset)
recent.grid(row = 1, column = 1, sticky = "WE")

ligprompt = tkinter.Label(main, text = "ligand chains ex: \"U+V\"")
ligprompt.configure(font = fontset)
ligprompt.grid(row = 2, column = 0, sticky = "W")
ligent = tkinter.Entry(main, bd = 5, textvariable = c)
ligent.configure(font = fontset)
ligent.grid(row = 2, column = 1, sticky = "WE")

rangeprompt = tkinter.Label(main, text = "file range inclusive ex: \"1:20\"")
rangeprompt.configure(font = fontset)
rangeprompt.grid(row = 3, column = 0, sticky = "W")
rangeent = tkinter.Entry(main, bd = 5, textvariable = d)
rangeent.configure(font = fontset)
rangeent.grid(row = 3, column = 1, sticky = "WE")

exclueprompt = tkinter.Label(main, text = "exclude files ex: \"5:10,2\"")
exclueprompt.configure(font = fontset)
exclueprompt.grid(row = 4, column = 0, sticky = "W")
exclueeent = tkinter.Entry(main, bd = 5, textvariable = e)
exclueeent.configure(font = fontset)
exclueeent.grid(row = 4, column = 1, sticky = "WE")

radprompt = tkinter.Label(main, text = "interaction radius in Angstroms")
radprompt.configure(font = fontset)
radprompt.grid(row = 5, column = 0, sticky = "W")
radent = tkinter.Entry(main, bd = 5, textvariable = f)
radent.configure(font = fontset)
radent.grid(row = 5, column = 1, sticky = "WE")

resprompt = tkinter.Label(main, text = "results directory (optional) ex: C:/users/user/.../folder")
resprompt.configure(font = fontset)
resprompt.grid(row = 6, column = 0, sticky = "W")
resent = tkinter.Entry(main, bd = 5, textvariable = g)
resent.configure(font = fontset)
resent.grid(row = 6, column = 1, sticky = "WE")
setres = tkinter.Button(main, text = "Choose", command = lambda: choosefile(0))
setres.grid(row = 6, column = 2, sticky = "W")

wkdirprompt = tkinter.Label(main, text = ".pdb files directory ex: C:/users/user/.../folder")
wkdirprompt.configure(font = fontset)
wkdirprompt.grid(row = 7, column = 0, sticky = "W")
wkdirent = tkinter.Entry(main, bd = 5, textvariable = h)
wkdirent.configure(font = fontset)
wkdirent.grid(row = 7, column = 1, sticky = "WE")
setwkdir = tkinter.Button(main, text = "Choose", command = lambda: choosefile(1))
setwkdir.grid(row = 7, column = 2, sticky = "W")

setprompt = tkinter.Label(main, text = "load settings from file ex: C:/users/user/.../folder")
setprompt.configure(font = fontset)
setprompt.grid(row = 8, column = 0, sticky = "W")
setent = tkinter.Entry(main, bd = 5, textvariable = i)
setent.configure(font = fontset)
setent.grid(row = 8, column = 1, sticky = "WE")
setbut = tkinter.Button(main, text = "Choose", command = iset)
setbut.grid(row = 8, column = 2, sticky = "W")
setload = tkinter.Button(main, text = "Load", command = load)
setload.grid(row = 8, column = 3)

savebut = tkinter.Checkbutton(main, text = "Do not save settings as new file", variable = save, onvalue = 1, offvalue = 0)
savebut.configure(font = fontset)
savebut.grid(row = 9, column = 1)

clear = tkinter.Button(main, text = "Clear", bd = 5, command = clr).grid(row = 10, column = 0, sticky = "E")
done = tkinter.Button(main, text = "Done", bd = 5, command = write).grid(row =10, column = 1, sticky = "W")

main.mainloop()


#post GUI processing
suffix = a.get()
chrec = b.get()
chlig = c.get()
intcut = float(f.get())

min = int(d.get().split(":")[0])
max = int(d.get().split(":")[1])
masterlist = list(range(min, max + 1))

store = e.get().split(",")
for i in range(0, len(store)):
	if re.match(".:.", store[i]):
		low = int(store[i].split(":")[0])
		high = int(store[i].split(":")[1])
		for j in range(low, high + 1):
			try:
				masterlist.remove(j)
			except ValueError:
				pass
	else:
		try:
			masterlist.remove(int(store[i]))
		except ValueError:
			pass

def surfaceCheck(iteration):
	cmd.get_sasa_relative("receptor." + str(iteration), var = "b", vis = 0, quiet = 1)

	cmd.select("check", selection = "chain " + chrec + " and b>0.7")
	modeloverlap = cmd.get_model("intrec." + str(iteration) + " in check")
	modelrec = cmd.get_model("intrec." + str(iteration))

	try:
		if(len(modeloverlap.atom) == len(modelrec.atom)):
			return True
		else:
			return False
	except ZeroDivisionError:
		return False

	cmd.remove("check")

#Load complexes and display as cartoon.
linew = "" #assign value as empty string before additions can be made
if (complete.get() == 1):
	
	
	for i in masterlist: #uses masterlist of all entries to be processed
		linew += str(i) + "    " #stores model number at left side then applys a tab
		fname = suffix + '.' + str(i) #stores as a string the filename to load into pymol
		cmd.load(h.get() + "/" + fname + '.pdb') #loads in fname model
		cmd.show("cartoon", fname) #displays model as cartoon
	
		#Choose receptor and ligand
		cmd.select("receptor." + str(i), selection = fname + " and chain " + chrec) #selects receptor by the name receptor.[modelnumber]
		cmd.select("ligand." + str(i), selection = fname + " and chain " + chlig) #selects ligand by the name ligand.[modelnumber]

		#Get interface
		cmd.select("intrec." + str(i), selection = "byres receptor." + str(i) + " within " + str(intcut) + " of ligand." + str(i)) #selects receptor interface by all residues that have atoms within [intcut]
		cmd.select("intlig." + str(i), selection = "byres ligand." + str(i) + " within " + str(intcut) + " of receptor." + str(i)) #selects ligand interface by all residues that have atoms within [intcut]
	
		go = surfaceCheck(i)
	
		#Write residue list to linew store
		if go:
			for sel in ["intrec." + str(i), "intlig." + str(i)]:
				reslist = ""
				model = cmd.get_model(sel)
				for obj in range(len(model.atom)):
					if model.atom[obj-1].resi != model.atom[obj].resi: #is there a posibility of missing any residues here???
						reslist += "%s%s%s," % (model.atom[obj].resn, model.atom[obj].resi, model.atom[obj].chain)
				linew += reslist + "    "
			linew += "\n"
		else:
			linew = linew + "invalid\n"

#OUTPUTFILE DETAILS
#Write interface residues to file for each complex
	f = open(g.get() + "/" + now + "_output.txt", 'w')
	header = "Format: Complex*Receptor interface*Ligand interface\n"
	f.write(header)
	f.write(linew)
	f.close()
	print("Program Complete")
	os.startfile(g.get() + "/" + now + "_output.txt")

else:
	print("Program Terminated")