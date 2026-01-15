# (c) 2026 IPK. All Rights Reserved
# Code written by: Amanda Souza Camara (camara@ipk-gatersleben.de)
# Based on Maksim Imakaev's scripts (from mirny lab)
from __future__ import print_function
from __future__ import division
import sys
import numpy as np 
from openmmlib_1 import Simulation
import polymerutils
import os
import simtk.unit as units
import pickle

nm = units.meter * 1e-9
fs = units.second * 1e-15
ps = units.second * 1e-12

N = int(sys.argv[1])   # number of monomers
# creating a starting conformation
oneVolume = 16000 #178479 #1331 #the volume occupied by one nucleosome, in cubic nm
totalVolume = N*oneVolume
# spherical confinement
radius = (totalVolume*3/(4*np.pi))**(1/3)

data = polymerutils.create_random_walk(1, N) # Creates an extended "random walk" conformation of length N
block=0
folder = str(sys.argv[3])
if not os.path.exists(folder):
 os.mkdir(folder)
 print("Directory", folder, " created.")
else:
 print("Directory ", folder, "already exists")

pi=np.pi

# Reading the distribution of centromeric and spacing regions
regionsFile = str(sys.argv[2])
file = open(regionsFile)
lines = file.readlines()
first = []
last = []
amount = []
for l in lines:
    first.append(int(l.split(' ')[0]))
    last.append(int(l.split(' ')[1]))
    amount.append(int(l.split(' ')[2]))
file.close()

centromeres = []
anchorParticles = []
cengroups = []
print("Creating random centromeric units according to each region.")
for i in range(len(first)):
    size =  last[i]-first[i]+1
    ncen = int(np.floor(amount[i]*size/100))
    if amount[i]>1: cunits = np.sort(first[i]+np.random.choice(np.arange(0,size),replace=False,size=(ncen))).tolist()
    if amount[i]==1: cunits = [np.int_(first[i]+size*np.random.random(amount[i]))]
    if amount[i]!=0:
        cengroups.append(cunits)
        centromeres = centromeres + cunits
        anchorParticles = anchorParticles + np.arange(cunits[0],cunits[-1]+1).tolist()
ncengroups=len(cengroups)
print (len(centromeres))
ncentromeres = len(centromeres)
nAnchorP = len(anchorParticles)
censtate = np.zeros(N)
for i in centromeres: censtate[i]=1


print("this simulation accounts for %d nucleosomes, %d centromeric groups, %d centromeric and %d anchor nucleosomes." % (N,ncengroups,ncentromeres,nAnchorP))

for i in range(0,totalBlocks,10000): saveBlock[i]=True
# Periodically tethering the centromeric nucleosomes to the center of mass of the centromeric unit.
for i in range(50000,totalBlocks,5): tetherBlock[i]=True
   
# saving special monomers final indexes to file
centromeres = np.array(centromeres)
pickle.dump(centromeres, open(os.path.join(folder, "centromeres.dat"),'wb'),protocol=2)
anchorParticles = np.array(anchorParticles)
pickle.dump(anchorParticles,open(os.path.join(folder,"anchors.dat"),'wb'),protocol=2)
#cengroups = np.array(cengroups)
pickle.dump(cengroups, open(os.path.join(folder, "cengroups.dat"),'wb'),protocol=2)

Epot = []
time = []
Ehand = []

print("Extrusion block steps are: %d" % (extrusionSteps))
blockRounds = firstBlocks+extrusionSteps++finalBlocks
print("The entire process will be over after %d blocks" % blockRounds)
totalParts = blockRounds//partialBlocks


# run the loop extrusion simulation in parts so it won't slow down with too many initialized bonds.
for part in range(totalParts):
 print("Starting part number %d" % (part))

 a = Simulation(thermostat=0.01)
 a.setup(platform="cuda", integrator="variableLangevin", errorTol=0.001, precision="mixed")
 a.load(data)  # loads a polymer, puts a center of mass at zero
 a.saveFolder(folder)

 a.addHarmonicPolymerBonds(wiggleDist=1, bondLength=10) # the same proportion used in the Maksim's loopExtrusion script
 
 #for c in cengroups:
 if block>=50000:
  k = float(sys.argv[4])
  a.tetherParticles(particles=centromeres,k=k, positions='current')
  print('Tethering %d particles.' % a.forceDict["Tethering Force"].getNumParticles())

 a.addPolynomialRepulsiveForce(trunc=5, radiusMult=10.5)
 a.addGrosbergStiffness(k=1.5)
 a.addSphericalConfinement(r=radius,k=5.)
# K is more or less arbitrary, k=4 corresponds to presistence length of 4,
# k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff

 numBonds = a.forceDict["HarmonicBondForce"].getNumBonds()
 print("There are %d polymer bonds." % numBonds)

# -----------Running a simulation ---------
 a.energyMinimization(stepsPerIteration=100)
 a.step = block
 if block==0: a.save()
 for b in range(partialBlocks):  # Do #blockRounds blocks
  
  if tetherBlock[block]:
   state = a.context.getState(getPositions=True)
   currentcoords = state.getPositions(asNumpy=True)
   comgroups=[]
   for c in cengroups: comgroups.append(np.sum(currentcoords[c],axis=0)/len(c)/nm)
   p=-1
   for g in range(ncengroups):
    for c in range(len(cengroups[g])):
     p=p+1
     a.forceDict["Tethering Force"].setParticleParameters(p,cengroups[g][c],comgroups[g])
     a.forceDict["Tethering Force"].updateParametersInContext(a.context)
  a.doBlock(100)  # Of 500 timesteps each
  if saveBlock[block]: a.save()

  Epot.append(a.state.getPotentialEnergy() / a.N / a.kT)
  time.append(a.state.getTime() / ps)
  
  data = a.getData()
  block = a.step

a.save()
outfile = open('%s/Epot.txt' % folder ,'w')
for i in range (len(Epot)):
 outfile.write("%5.2f\n" % Epot[i])
outfile.close()
