import pandas as pd
import numpy as np

df = pd.read_csv('array_sizes_and_spacings_20_species.csv')
species=df['genome'].unique().tolist()

for s in species:
    print(s)
    m=df.loc[df['genome'] == s][['array_size','following_spacing_size']]
    narrays=len(m)
    sums=[]
    for i in range(narrays-1):
        for j in range(i+1,narrays):
            if not m.iloc[i:j,:].isnull().any().any(): sums.append([i,j,np.sum(np.sum(m.iloc[i:j,:],axis=0),axis=0),j-i])
    sumsdf = pd.DataFrame(sums, columns = ['start', 'end', 'sum','narrays'])
    print(sumsdf)

    # take the three regions whose sum of array and space sizes are the closest to 15Mb
    sumsdf_15Mb = sumsdf.loc[sumsdf['sum']<15000000]
    max3_sum=sumsdf_15Mb.loc[sumsdf_15Mb['sum'].nlargest(3).index]
    print(max3_sum)
    for nmax in range(3):
        filename=s+'_regions'+str(nmax+1)+'.txt' # save each pf the three options in separate files
        outfile = open(filename,'w')
        max_sum=max3_sum.iloc[nmax,:]
        first_spacing=15000000-int(max_sum['sum'])
        outfile.write("%d %d %d \n" % (0,15000000-int(max_sum['sum']),0))
        start=first_spacing
        for i in range(int(max_sum['start']),int(max_sum['end'])):
            array_size=int(m.iloc[i,0])
            following_spacing_size=int(m.iloc[i,1])
            outfile.write("%d %d %d \n" % (start//200,(start+array_size-1)//200,60)) # divide by 200 
            outfile.write("%d %d %d \n" % ((start+array_size)//200,(start+array_size+following_spacing_size-1)//200,0))
            start=start+array_size+following_spacing_size
        outfile.close()
