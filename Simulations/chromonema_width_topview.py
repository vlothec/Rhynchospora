import numpy as np
import sys
import joblib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


filename = str(sys.argv[1])
file = open(filename)
fileslist = file.readlines()
file.close()
datafiles = []
for f in fileslist: datafiles.append(f.rstrip('\n'))
speciesdict = {}
for d in datafiles:
    species=d.rsplit('/')[1]
    speciesdict[species]=[]

filename = str(sys.argv[2])
file = open(filename)
fileslist = file.readlines()
file.close()
regionsfiles = []
for f in fileslist: regionsfiles.append(f.rstrip('\n'))

#colors=['#ef476f','#ef476f','#ef476f','#f78c6b','#f78c6b','#f78c6b','#ffd166','#ffd166','#ffd166','#06d6a0','#06d6a0','#06d6a0','#118ab2','#118ab2','#118ab2','#073b4c','#073b4c','#073b4c']
colors=['#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c','#ef476f','#f78c6b','#ffd166','#06d6a0','#118ab2','#073b4c']


#datafiles=['species_results/breviuscula/1/block50000.dat','species_results/pubera/1/block50000.dat','species_results/rugosa/1/block50000.dat']
#regionsfiles=['regions_modeller/regions_breviuscula1.txt','regions_modeller/regions_pubera1.txt','regions_modeller/regions_rugosa1.txt']
#colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


nfiles=len(datafiles)
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
    for l in lines:
        tag=int(l.strip('\n').split()[2])
        #print(l.strip('\n').split()[2])
        if tag==0:
            regions.append(range(int(l.strip('\n').split()[0]),int(l.strip('\n').split()[1])))
            type.append(tag)
        if tag!=0:
            regions.append(range(int(l.strip('\n').split()[0]),int(l.strip('\n').split()[1])))
            type.append(tag)


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
            #for j in range(len(data_c)): speciesdict[species].append(np.sqrt(np.sum(np.cross(data_c[j],a)**2))/moda)
            for j in range(len(data_c)): speciesdict[species].append(abs(np.dot(data_c[j],ref)/modref))
            #width.append(np.mean(proj))
    #print(len(proj))

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
#initial_guess=[120,80,80,
#               100,100,80]
for s in range(len(speciesdict)):
#for s in [0]:
    species = list(speciesdict.keys())[s]
    n,bins=np.histogram(speciesdict[species],bins=50)
    bin_centres = (bins[:-1] + bins[1:])/2.

    ydata = n
    xdata = bin_centres
    p0 = [max(ydata)-min(ydata), initial_guess[s],1,min(ydata)] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
    print(species, popt[1])

    print(xdata[list(sigmoid(xdata, *popt)).index(min(sigmoid(xdata, *popt), key=lambda x:abs(x-1000)))])
    #print(xdata,sigmoid(xdata, *popt))
    plt.plot(xdata, sigmoid(xdata, *popt), color=colors[s])
    plt.scatter(bin_centres, n, marker='o', s=40, alpha=1,label=species,color=colors[s])
    #plt.axvline(x = np.mean(speciesdict[species]), color=colors[s])
    #print(species,np.mean(speciesdict[species]))
    plt.title(species)
    plt.show()

#plt.xscale('log')
plt.xlabel('Radius (nm)')
plt.ylabel('Number of nucleosomes (#)')
plt.legend()
plt.gcf().set_size_inches(10, 5)
#plt.show()
plt.savefig('chromonema_width_topview.pdf',dpi=300, bbox_inches='tight')






#print(np.mean(width))





#
#    data_c=data[cens[i]]-com
#    rg.append(np.sqrt(np.sum(data_c**2)/lcen))
#    print('c',rg[i])
#print(np.mean(rg),np.std(rg))
    # for control
#    controlsel=[j+200 for j in sels[i]]
#    lcontrolsel=len(controlsel)
#    com = np.sum(data[controlsel],axis=0)/lcontrolsel
#    data_c=data[controlsel]-com
#    rg = np.sqrt(np.sum(data_c**2)/lcontrolsel)
#    print(rg)
