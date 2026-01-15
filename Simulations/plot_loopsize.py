import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.signal import savgol_filter

filename = str(sys.argv[1]) # a file containing a list of paths for the loopsize.dat files output from simulations.
samplename = filename.split('/')[1]
file = open(filename)
lines = file.readlines()
file.close()

loopsize=[]
for l in lines: loopsize.append(float(l.rstrip('\n')))
print(samplename," Mean loop size in the end: %.2f" % np.mean(loopsize[-50:]))
x=np.arange(len(loopsize))

loopsize_smooth = savgol_filter(loopsize,window_length=20000, polyorder=2)

plt.plot(loopsize, label='raw', alpha=0.5)
plt.plot(loopsize_smooth, label='smoothed', linewidth=2)
plt.legend()
plt.savefig(samplename+'_loopsize.png', dpi=300)
#plt.show()
plt.close()

derivative=np.gradient(loopsize_smooth,x)
derivative_smooth = savgol_filter(derivative,window_length=20000, polyorder=2)

for i, val in enumerate(derivative_smooth):
    if val <= 0: 
        print(samplename,' First negative derivative at',i,val)
        break

plt.plot(derivative,label='derivative of smoothed loopsize')
plt.plot(derivative_smooth,label='smoothed derivative')
plt.savefig(samplename+'_derivativeloopsize.png', dpi=300)
#plt.show()
plt.close()