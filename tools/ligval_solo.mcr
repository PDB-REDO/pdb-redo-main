# ligval_solo.mcr: a YASARA macro to validate ligands and their interactions.
# Version 1.00
#
# Minimum YASARA tier: YASARA model, YASARA dynamics is required to calculate the energy of formation.
#
# This script was created by Robbie P. Joosten (r.joosten@nki.nl, robbie_joosten@hotmail.com)
# Reference: If you publish results (directly or indirectly) obtained by using this macro, please cite YASARA and (any of)
# these publications:
# 1) Adria Cereto-Massague, Maria J. Ojeda, Robbie P. Joosten, Cristina Valls, Miquel Mulero, M. Josepa Salvado, 
#    Anna Arola-Arnal, Lluis Arola, Santiago Garcia-Vallvé, Gerard Pujadas: "The good, the bad and the dubious: VHELIBS, a 
#    validation helper for ligands and binding sites" J Cheminform. 5(1) p. 36 (2013)
# 2) Robbie P. Joosten, Fei Long, Garib N. Murshudov, Anastassis Perrakis: "The PDB_REDO server for macromolecular
#    structure model optimization" IUCrJ 1, p. 213-220 (2014)
#
# Changelog:
# Version 1.00
# - First attempt based on the new versus old ligand validation script.
#
# YASARA input variables:
# newpdb  The PDB file of the changed model
# resnum  The residue number of the ligand to be inspected
# chid    The chain ID of the ligand to be inspected
#
#Initialise
OnError Exit

#Check the YASARA version
if (Model)==0
  Print "YASARA version too low, you need at least YASARA Model"
  exit
  
Console off


#Load the files
LoadPDB (newpdb)

#Add hydrogens
AddHydAll

#Count the bumps
print "Start bumps new"
ListConRes Res (resnum) Mol (chid) Obj 1, Obj 1, Cutoff=-0.42, Subtract=HBoRadii, Exclude=5, Sort=No
print "End bumps new"

#Get the hydrogen bonds and their energy
print "Start hbond new"
ListHBoRes Res (resnum) Mol (chid) Obj 1, Obj 1
print "End hbond new"

#Get hydrophobic inteactions
print "Start hydpho new"
ListIntRes Res (resnum) Mol (chid) Obj 1, Obj 1, Type=Hydrophobic, Exclude=5, Occluded=No, Sort=No
print "End hydpho new"

#Get pi-pi interactions
print "Start pipi new"
ListIntRes Res (resnum) Mol (chid) Obj 1, Obj 1, Type=PiPi, Exclude=5, Occluded=No, Sort=No
print "End pipi new"

#Get cation-pi interactions
print "Start catpi new"
ListIntRes Res (resnum) Mol (chid) Obj 1, Obj 1, Type=CationPi, Exclude=5, Occluded=No, Sort=No
print "End catpi new"


#Check the YASARA again
if (Dynamics)==0
  Print "YASARA version too low to calculate ligand heat of formation, you need at least YASARA dynamics"
  exit
  
#Calculate formation energy
FormEnergyRes Res (resnum) Mol (chid) Obj 1

exit
