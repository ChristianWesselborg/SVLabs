#Delete previously loaded molecules just in case
delete all

#Change directory to current directory
cd c:/Users/WESSEC/Documents/Research/RPI GB/pdbfiles/

#Load molecule
load 6ira.pdb

#start from scratch 
hide everything 

as cartoon, 6IRA
util.cbc 6IRA

#SELECT TMD 
sel tmdA, 6IRA and chain A and ( resid 546:657 or resid 802:846)
sel tmdB, 6IRA and chain B and (resi 549:657 or resi 807:841)
sel tmdC, 6IRA and chain C and ( resid 546:657 or resid 802:846)
sel tmdD, 6IRA and chain D and (resi 549:657 or resi 807:841)
sel tmd, tmdA or tmdB or tmdC or tmdD
as surface, tmd
set transparency, 0.75

#Block residue selection
sel blockA, chain A and (resid 544+545+658+659+800+801)
sel blockB, chain B and (resid 547+548+658+659+806+807)
sel blockC, chain C and (resid 544+545+658+659+800+801)
sel blockD, chain D and (resid 547+548+658+659+806+807)

#deletions and representations 
as sphere, blockA
as sphere, blockB
as sphere, blockC
as sphere, blockD
sel blocked, blockA or blockB or blockC or blockD

color gray, blobked
set sphere_scale, 2.5
color gray, blocked

sel rest, not tmd
set bg_rgb, [ 1, 1, 1]

#Save files
#save 6IRA_notmd.pdb, rest #uncomment for use
