from math import *
from random import *
import numpy as np

ParticleN = len(np.loadtxt('ThreeLetterCode.dat', dtype="U", unpack=True))  # Determine the length of the IDP
print ("IDP Length:", ParticleN)

fi = open (str(ParticleN),"w")

spacex = 0.38

x=[]
y=[]
z=[]


x.append(10.0)
y.append(0.0)
z.append(0.0)

i = 1
while (i<ParticleN):
    count =0
    xi = uniform(0,0.02)
    #yi = uniform(0,0.01)
    yi = 0.0
    x.append (x[0]+spacex*i)
    y.append (y[0]+yi)
    z.append (0.0)
    i=i+1

for i in range (ParticleN):
    print (x[i],y[i],z[i])
    print (x[i],y[i],z[i], file = fi)


for i in range (ParticleN):
    for j in range (ParticleN):
        if (i != j):
            r = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))
            if (r < spacex-0.1):
                print ("------------------------------------------------")
                print (i, j, r)
                
                
# Input parameters for all the amino acids (force field)
ff_para = 'stats_module.dat'
aalist={}
with open(ff_para,'r') as fid:
	for i in fid:
		if i[0]!='#':
			tmp=i.rsplit()
			aalist[tmp[0]]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())

# This translates each amino acid type into a number, which will be used in HOOMD
# For example, GLY is with an ID of 10

aamass=[]
aacharge=[]
aasigma=[]
aahps=[]
print ("List of Amino Acids:", aakeys)

for i in aakeys:
	aamass.append(aalist[i][0])
	aacharge.append(aalist[i][1])
	aasigma.append(aalist[i][2])
	aahps.append(aalist[i][3])

# Now we can translate the entire sequence into a number code according to the order in 'aakeys'

fi_param = open ("chain_param.dat","w")
with open('ThreeLetterCode.dat','r') as fid:
    for i in fid:
        iname=i.rsplit()[0]
        print(iname, aakeys.index(iname), aalist[iname][0], aalist[iname][1], aalist[iname][2], aalist[iname][3], file = fi_param)
        
fid.close()
fi_param.close()
















