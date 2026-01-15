# (c) 2026 IPK. All Rights Reserved
# Code written by: Amanda Souza Camara (camara@ipk-gatersleben.de)
# Based on Maksim Imakaev's scripts (from mirny lab)
from __future__ import print_function
from __future__ import division
import sys
import numpy as np 
from openmmlib_1 import Simulation
import os
import simtk.unit as units
import bondsPredictor11
import pickle
import joblib

nm = units.meter * 1e-9
fs = units.second * 1e-15
ps = units.second * 1e-12

N = int(sys.argv[1])   # number of monomers
# parameters for the SMC bonds.
nLEFI = int(sys.argv[3]) # this is related to separation of condensins I
nLEFII = int(sys.argv[4])
nLEFs = nLEFI + nLEFII
# Parameter for binding two nucleosomes as if they were held by a loop extruder.  The same proportion used in the loopExtrusion script from Maksim.
wiggleDist = 2 # Average displacement from the equilibrium bond distance
bondLength = 5 # The length of the SMC bonds
# bond dictionaries parameters
extrusionSteps = 50000
consideredSteps = 50000
firstBlocks = 0 # how many blocks to wait before loop extrusion begins
waitBlocks = 0 # how many blocks to wait after extrusion of condensins II to add condensins I
partialBlocks = 50 # how many blocks per initialized simulation
prophaseBlocks = 0 # how many blocks to run after all loops of condensin I are extruded
finalBlocks = 0
totalBlocks = firstBlocks + extrusionSteps + finalBlocks
totalParts = int(totalBlocks/partialBlocks)
saveBlock = [False]*(totalBlocks+1) 
tetherBlock = [False]*(totalBlocks+1)

# uploads a conformation from file
blockofile = sys.argv[2]
backData = joblib.load(blockofile)
ifpos = blockofile.rfind('/')+1
infolder = blockofile[:ifpos]
print('input folder', infolder)
data = np.array(backData['data'])
block=0
folder = str(sys.argv[7])
if not os.path.exists(folder):
 os.mkdir(folder)
 print("Directory", folder, " created.")
else:
 print("Directory ", folder, "already exists")

pi=np.pi

# Reading data on centromeres created in the previous simulation without loop extrusion
centromeres = joblib.load(infolder+'centromeres.dat')
anchors = joblib.load(infolder+'anchors.dat')
cengroups = joblib.load(infolder+'cengroups.dat')
print (len(centromeres))
ncentromeres = len(centromeres)
nanchors = len(anchors)
ncengroups = len(cengroups)

print("this simulation accounts for %d nucleosomes, being %d centromeric, %d anchor nucleosomes and %d condensins I and %d condensis II." % (N,ncentromeres,nanchors,nLEFI,nLEFII))

for i in range(0,totalBlocks,10000): saveBlock[i]=True

# Reduce tethering to less steps
for i in range(0,5000,5): tetherBlock[i]=True
for i in range(5000,10000,10): tetherBlock[i]=True
for i in range(10000,totalBlocks,50):tetherBlock[i]=True

# list of lifetimes
lifetime = []
for i in range(nLEFII): lifetime.append(int(sys.argv[6]))
for i in range(nLEFI): lifetime.append(int(sys.argv[5]))

#list of first binding time
birth = []
for i in range(nLEFII): birth.append(1)
for i in range(nLEFI): birth.append(100000000)


# list of LEFs being anchored by centromeres
centromereAnchor=[]
for i in range(nLEFII): centromereAnchor.append(1)
for i in range(nLEFI): centromereAnchor.append(0)

# predicting loop dynamics with fortran   
passElements = consideredSteps*nLEFs*2
centromeres = np.array(centromeres)+1
anchors = np.array(anchors)+1
bondDict = bondsPredictor11.extruder(N,passElements,extrusionSteps,lifetime,birth,centromereAnchor,anchors,nLEFs,nanchors)
centromeres = centromeres - 1
anchors = anchors - 1
extrusionSteps = consideredSteps
print(bondDict, bondDict.shape)
loopBondDict = np.zeros((extrusionSteps,nLEFs,5))
loopBondDict[:,:,0:2] = bondDict.reshape(extrusionSteps,nLEFs,2) - 1
loopBondDict = loopBondDict.astype(int)
print(loopBondDict.shape)
print("average loop size is %d." % np.average([loopBondDict[extrusionSteps-1,:,1]-loopBondDict[extrusionSteps-1,:,0]]))
outfile = open("%s/loopsize.dat" % folder, 'w')
for i in range(extrusionSteps):
 outfile.write("%5.2f\n" % np.average(loopBondDict[i,:nLEFII,1]-loopBondDict[i,:nLEFII,0]))
outfile.close()

# saving special monomers final indexes to file
centromeres = np.array(centromeres)
pickle.dump(centromeres, open(os.path.join(folder, "centromeres.dat"),'wb'),protocol=2)
anchors = np.array(anchors)
pickle.dump(anchors,open(os.path.join(folder,"anchors.dat"),'wb'),protocol=2)
condensinI = loopBondDict[extrusionSteps-1,nLEFII:,0:2]
pickle.dump(condensinI, open(os.path.join(folder, "condensinI.dat"),'wb'),protocol=2)
condensinII = loopBondDict[extrusionSteps-1,:nLEFII,0:2]
pickle.dump(condensinII, open(os.path.join(folder, "condensinII.dat"),'wb'),protocol=2)


# do the loop extrusion => activates and inactivates bonds changing the value of k.
class smcMilker:
 def step(self, context, lBD, part, b, k, smcBondLength):
  print("milking",part,b)
  blockNumber = part*partialBlocks + b
  cIIfile = "block" + str(blockNumber) + "condensinII.dat"
  if (saveBlock[blockNumber]):
   pickle.dump(lBD[blockNumber,:nLEFII,0:2], open(os.path.join(folder, cIIfile),'wb'),protocol=2)
  for i in range(extrusionSteps):
   if lBD[i,0,2]==part and lBD[i,0,3]==b :
    for h in range (nLEFs):
     if lBD[i,h,0]!=-1: a.forceDict["HarmonicBondForce"].setBondParameters(int(lBD[i,h,4]), int(lBD[i,h,0]), int(lBD[i,h,1]), smcBondLength, k)
   elif lBD[i,0,2]==part and lBD[i,0,3]!=b:
    for h in range (nLEFs):
     if lBD[i,h,0]!=-1:a.forceDict["HarmonicBondForce"].setBondParameters(int(lBD[i,h,4]), int(lBD[i,h,0]), int(lBD[i,h,1]), smcBondLength, 0)
  a.forceDict["HarmonicBondForce"].updateParametersInContext(context)

 def initialize(self, lBD,part,k,smcBondLength):
  print('setting up SMC bonds for this part',part )
  for b in range (extrusionSteps):
   if lBD[b,0,2]==part:
    for h in range (nLEFs):
     if lBD[b,h,0]!=-1: lBD[b,h,4] = a.forceDict["HarmonicBondForce"].addBond(int(lBD[b,h,0]),int(lBD[b,h,1]),smcBondLength,0)
  numBonds = a.forceDict["HarmonicBondForce"].getNumBonds()
  print("This simulation accounts now for %d bonds" % numBonds)

milker = smcMilker()

Epot = []
time = []
Ehand = []

print("Extrusion block steps are: %d" % (extrusionSteps))
blockRounds = firstBlocks+extrusionSteps++finalBlocks
print("The entire process will be over after %d blocks" % blockRounds)
totalParts = blockRounds//partialBlocks
for i in range(extrusionSteps):
 currentPart = (firstBlocks+i)//partialBlocks
 currentBlock = int((firstBlocks+i)%partialBlocks)
 loopBondDict[i,:,2] = currentPart
 loopBondDict[i,:,3] = currentBlock


# run the loop extrusion simulation in parts so it won't slow down with too many initialized bonds.
for part in range(totalParts):
 print("Starting part number %d" % (part))

 a = Simulation(thermostat=0.01)
 a.setup(platform="cuda", integrator="variableLangevin", errorTol=0.001, precision="mixed")
 a.load(data)  # loads a polymer, puts a center of mass at zero
 a.saveFolder(folder)

 a.addHarmonicPolymerBonds(wiggleDist=1, bondLength=10) # the same proportion used in the Maksim's loopExtrusion script
 a.addPolynomialRepulsiveForce(trunc=5, radiusMult=10.5)
 a.addGrosbergStiffness(k=0.5)

 a.tetherParticles(particles=centromeres,k=0.01, positions='current')
 print('Tethering %d particles.' % a.forceDict["Tethering Force"].getNumParticles())

 k = a.kbondScalingFactor / (wiggleDist ** 2) # harmonic force constant
 smcBondLength = bondLength * a.length_scale

 numBonds = a.forceDict["HarmonicBondForce"].getNumBonds()
 print("There are %d polymer bonds." % numBonds)

 milker.initialize(loopBondDict,part,k,smcBondLength)
 numBonds = a.forceDict["HarmonicBondForce"].getNumBonds()
 print(numBonds)

# -----------Running a simulation ---------
 a.energyMinimization(stepsPerIteration=100)
 a.step = block
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

  milker.step(a.context,loopBondDict,part,b,k,smcBondLength)
 
  a.doBlock(100)  # Of 100 timesteps each
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
