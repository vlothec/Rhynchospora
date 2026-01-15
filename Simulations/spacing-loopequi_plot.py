import numpy as np
import sys
import matplotlib.pyplot as plt

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


filename = str(sys.argv[1]) # regions files (output of find arrays)
file = open(filename)
fileslist = file.readlines()
file.close()
regionsfiles = []
for f in fileslist: regionsfiles.append(f.rstrip('\n'))

filename = str(sys.argv[2]) # loop size equilibration time file
file = open(filename)
fileslist = file.readlines()
file.close()
loopequi = {}
for f in fileslist:
    loopequi[f.rstrip('\n').split(' ')[0]]=int(f.rstrip('\n').split(' ')[1])*100/1000000
print(loopequi)

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


nfiles=len(regionsfiles)
spacing_mean = {}
spacing_min = []
spacing_max = []
looplength_mean =[]
looplength_std = []
samples = []
for f in range(nfiles):

    regionsfile=regionsfiles[f]
    sample = regionsfile.split('/')[1][:-12]+'rep'+regionsfile.split('/')[1][-5:-4]
    samples.append(sample)
    if not sample in loopequi.keys(): loopequi[sample]= 50000*100/1000000
    file = open(regionsfile)
    lines = file.readlines()
    file.close()

    nregions=len(lines)
    regions = []
    type = []
    spacings = []
    for l in lines:
        tag=int(l.strip('\n').split()[2])
        if tag==0:
            regionmonomers=range(int(l.strip('\n').split()[0]),int(l.strip('\n').split()[1]))
            spacings.append(len(regionmonomers)*200/1000)
    spacing_mean[sample]=np.mean(spacings)
    spacing_min.append(np.min(spacings))
    spacing_max.append(np.max(spacings))
    print(samples[-1],spacing_mean[samples[-1]])

for i in range(len(samples)):
    s=samples[i][:-5]
    if s!='Rhync_filiformis_5519G.hap1.chr' and s!='R_barbata_hap1_chrs':
        plt.scatter(spacing_mean[samples[i]],loopequi[samples[i]],c=colors[s],label=s,marker=markers[s])
plt.xlabel('Spacing (kb)')
plt.ylabel('Loop length equilibration time (1 million simulation steps)')
plt.savefig('LoopEquilibration_Spacing.pdf', dpi=300)
plt.show()
