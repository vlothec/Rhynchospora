import numpy as np
import sys
import joblib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


filename = str(sys.argv[1]) # block files
file = open(filename)
fileslist = file.readlines()
file.close()
datafiles = []
for f in fileslist: datafiles.append(f.rstrip('\n'))
speciesdict = {}
speciesdict_centerview = {}
for d in datafiles:
    species=d.rsplit('/')[1]
    speciesdict[species]=[]
    speciesdict_centerview[species]=[]

filename = str(sys.argv[2]) # regions files
file = open(filename)
fileslist = file.readlines()
file.close()
regionsfiles = []
for f in fileslist: regionsfiles.append(f.rstrip('\n'))

filename = str(sys.argv[3]) # loop lengths files
file = open(filename)
fileslist = file.readlines()
file.close()
loopsfiles = []
for f in fileslist: loopsfiles.append(f.rstrip('\n'))

#colors=['#ef476f','#ef476f','#ef476f','#f78c6b','#f78c6b','#f78c6b','#ffd166','#ffd166','#ffd166','#06d6a0','#06d6a0','#06d6a0','#118ab2','#118ab2','#118ab2','#073b4c','#073b4c','#073b4c']
#colors=['#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c']
colors={'R_barbata_hap1_chrs':'#7470AE',
        'R_holoschoenoides_hap1_chr':'#6FB2E4',
        'R_pubera_ref_2n10_v2.chr':'#B0D767',
        'Ralba.chr':'#C66526',
        'Rbreviuscula.hap1.chr':'#7DC0A6',
        'Rcephalotes.hap2.chr':'#ED936B',
        'Rcolorata.hap1.chr':'#75A43A',
        'Rhync_austrobrasiliensis_6344D.hic.hap1.chr':'#D43E88',
        'Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr':'#4B9C7A',
        'Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr':'#3070AC',
        'Rhync_filiformis_5519G.hap1.chr':'#DA8EC0',
        'Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr':'#469C76',
        'Rhync_radicans.asm.bp.p_ctg.FINAL.chr':'#F9DA56',
        'Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr':'#E0C59A',
        'Rhync_riparia_5519C.hap1.chr':'#919FC7',
        'Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr':'#EEE461',
        'Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr':'#CA6627',
        'Rrugosa.chr':'#DCA237',
        'Rtenerrima.chr':'#C17DA5',
        'Rtenuis.hap1.chr':'#B3B3CA'}
markers={'R_barbata_hap1_chrs':'o',
        'R_holoschoenoides_hap1_chr':'^',
        'R_pubera_ref_2n10_v2.chr':'D',
        'Ralba.chr':'o',
        'Rbreviuscula.hap1.chr':',',
        'Rcephalotes.hap2.chr':'D',
        'Rcolorata.hap1.chr':'D',
        'Rhync_austrobrasiliensis_6344D.hic.hap1.chr':'^',
        'Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr':'^',
        'Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr':'D',
        'Rhync_filiformis_5519G.hap1.chr':'o',
        'Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr':'o',
        'Rhync_radicans.asm.bp.p_ctg.FINAL.chr':'^',
        'Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr':',',
        'Rhync_riparia_5519C.hap1.chr':',',
        'Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr':'o',
        'Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr':'^',
        'Rrugosa.chr':',',
        'Rtenerrima.chr':',',
        'Rtenuis.hap1.chr':'D'}


#datafiles=['species_results/breviuscula/1/block50000.dat','species_results/pubera/1/block50000.dat','species_results/rugosa/1/block50000.dat']
#regionsfiles=['regions_modeller/regions_breviuscula1.txt','regions_modeller/regions_pubera1.txt','regions_modeller/regions_rugosa1.txt']
#colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


nfiles=len(datafiles)
spacing_mean = []
spacing_min = []
spacing_max = []
looplength_mean =[]
looplength_std = []
for f in range(nfiles):
    datafile = datafiles[f]
    print(datafile)
    block=joblib.load(datafile)
    data = np.array(block['data'])
    #data = data - np.minimum(np.min(data, axis=0), np.zeros(3, float) - 100)[None, :]
    N=data.shape[0]
    d=data.shape[1]
    #print(N,d)

    species = datafile.rsplit('/')[1]

    regionsfile=regionsfiles[f]
    print(regionsfile)
    file = open(regionsfile)
    lines = file.readlines()
    file.close()

    nregions=len(lines)
    regions = []
    type = []
    spacings = []
    for l in lines:
        tag=int(l.strip('\n').split()[2])
        #print(l.strip('\n').split()[2])
        if tag==0:
            regionmonomers=range(int(l.strip('\n').split()[0]),int(l.strip('\n').split()[1]))
            regions.append(regionmonomers)
            type.append(tag)
            spacings.append(len(regionmonomers)*200/1000)
        if tag!=0:
            regions.append(range(int(l.strip('\n').split()[0]),int(l.strip('\n').split()[1])))
            type.append(tag)
    spacing_mean.append(np.mean(spacings))
    spacing_min.append(np.min(spacings))
    spacing_max.append(np.max(spacings))

    loopsfile=loopsfiles[f]
    print(loopsfile)
    file = open(loopsfile)
    lines = file.readlines()
    file.close()
    ll=[]
    for l in lines:
        ll.append(float(l.strip('\n').split()[0])*200/1000)
    #print(ll)
    looplength_mean.append(np.mean(ll[40000:]))
    looplength_std.append(np.std(ll[40000:]))

    #print(regions)
    #print(nregions)

    #width=[]
    #proj=[]
    for i in range(1,nregions-1):
        if type[i]==0:
            lreg0=len(regions[i-1])
            lreg1=len(regions[i])
            lreg2=len(regions[i+1])
            com0 = np.sum(data[regions[i-1]],axis=0)/lreg0
            com1 = np.sum(data[regions[i]],axis=0)/lreg1
            com2 = np.sum(data[regions[i+1]],axis=0)/lreg2
            data_c=data[regions[i]]-com1
            #print(i)
            #print(len(data_c))
            a = com2-com0
            moda = np.sqrt(np.sum(a**2))
            ref = np.cross(data_c[0],a)
            modref = np.sqrt(np.sum(ref**2))
            #proj=[]
            for j in range(len(data_c)): speciesdict_centerview[species].append(np.sqrt(np.sum(np.cross(data_c[j],a)**2))/moda)
            for j in range(len(data_c)): speciesdict[species].append(abs(np.dot(data_c[j],ref)/modref))
            #width.append(np.mean(proj))
    #print(len(proj))
sigthreshold1000=[]
siginflection=[]
meanradius=[]
initial_guess=[120,150,120,
               120,120,100,
               120,100,100,
               120,120,100,
               120,100,100,
               120,80,80,
               120,120,120,
               130,120,120,
               120,120,120,
               120,120,120,
               120,120,120,
               120,100,100,
               100,120,80,
               100,100,80,
               80,120,100,
               100,120,120,
               130,120,120,
               120,100,100,
               80,100,120,
               100,120,100]
for s in range(len(speciesdict)):
#for s in [0]:
    species = list(speciesdict.keys())[s]
    n,bins=np.histogram(speciesdict[species],bins=50)
    bin_centres = (bins[:-1] + bins[1:])/2.

    ydata = n
    xdata = bin_centres
    p0 = [max(ydata)-min(ydata), initial_guess[s],1,min(ydata)] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
    #print(species, popt[1])
    sigthreshold1000.append(popt[1])

    sig_inflection=xdata[list(sigmoid(xdata, *popt)).index(min(sigmoid(xdata, *popt), key=lambda x:abs(x-1000)))]
    #print(sig_inflection)
    siginflection.append(sig_inflection)

    meanradius.append(np.mean(speciesdict_centerview[species]))
    print(species,'mean radius ', meanradius[s],'mean spacing ', spacing_mean[s])

    #print(xdata,sigmoid(xdata, *popt))


#print(len(siginflection),len(sigthreshold1000),len(spacing_mean),len(spacing_max),len(spacing_min),len(looplength_mean),len(looplength_std),len(meanradius))
species=list(speciesdict.keys())

all_results = pd.DataFrame(
    {'sample':species,
     'spacing':spacing_mean,
     'loop_length':looplength_mean,
     'centerviwe_mean_radius':meanradius,
     'topview_inflection':siginflection,
     'topview_threshold':sigthreshold1000
     })
all_results.to_csv('simulation_results.csv')

for i in range(len(species)):
    s=species[i][:-5]
    if s!='Rhync_filiformis_5519G.hap1.chr' and s!='R_barbata_hap1_chrs':
        plt.scatter(spacing_mean[i],looplength_mean[i],c=colors[s],label=s,marker=markers[s])
plt.xlabel('Spacing (kb)')
plt.ylabel('Loop length (kb)')
plt.savefig('LoopLength_Spacing.pdf', dpi=300)
plt.show()

for i in range(len(species)):
    s=species[i][:-5]
    if s!='Rhync_filiformis_5519G.hap1.chr' and s!='R_barbata_hap1_chrs':
        plt.scatter(spacing_mean[i],meanradius[i],c=colors[s],label=s,marker=markers[s])
plt.xlabel('Spacing (kb)')
plt.ylabel('Mean Radius (nm)')
plt.savefig('MeanRadius_Spacing.pdf', dpi=300)
plt.show()

for i in range(len(species)):
    s=species[i][:-5]
    if s!='Rhync_filiformis_5519G.hap1.chr' and s!='R_barbata_hap1_chrs':
        plt.scatter(spacing_mean[i],siginflection[i],c=colors[s],label=s,marker=markers[s])
plt.xlabel('Spacing (kb)')
plt.ylabel('Thickness top view sigmoid inflection (nm)')
plt.savefig('ThicknessInf_Spacing.pdf', dpi=300)
plt.show()

for i in range(len(species)):
    s=species[i][:-5]
    if s!='Rhync_filiformis_5519G.hap1.chr' and s!='R_barbata_hap1_chrs':
        plt.scatter(spacing_mean[i],sigthreshold1000[i],c=colors[s],label=s,marker=markers[s])
plt.xlabel('Spacing (kb)')
plt.ylabel('Thickness top view sigmoid threshold (nm)')
plt.savefig('ThicknessThr_Spacing.pdf', dpi=300)
plt.show()
