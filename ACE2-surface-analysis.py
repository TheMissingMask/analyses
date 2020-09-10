#!/usr/bin/python

import numpy as np
import MDAnalysis
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import seaborn as sns

def runningMean(x,N):
    cumsum=np.cumsum(np.insert(x,0,0))
    return (cumsum[N:]-cumsum[:-N])/float(N)

u=MDAnalysis.Universe('centred.pdb','centred.xtc')
prot=u.select_atoms('protein')
sol=u.select_atoms('resname WN or resname CL or resname NA')
        
saDict={}
for ts in u.trajectory:
    saDict[ts.time]=[]
    for res in prot.residues:
        minD=np.min(cdist(res.atoms.positions,sol.atoms.positions))
        if minD<=5.0:
            saDict[ts.time].append(1)
        else:
            saDict[ts.time].append(0)
            
resList=[x+1 for x in range(len(saDict[0]))]
        
saList=[0]*len(resList)
for i in range(len(resList)):
    for k in saDict.keys():
        saList[i]+=saDict[k][i]
        
saList=[x/len(saDict.keys()) for x in saList]

sns.set_style('white')

resListWT=runningMean(resList,20)
saListWT=runningMean(saList,20)
print('done wt\n')

u=MDAnalysis.Universe('No90Glycan/MD/md.tpr','No90Glycan/MD/md.xtc')
prot=u.select_atoms('protein')
sol=u.select_atoms('resname WN or resname CL or resname NA')

saDict={}
for ts in u.trajectory:
    saDict[ts.time]=[]
    for res in prot.residues:
        minD=np.min(cdist(res.atoms.positions,sol.atoms.positions))
        if minD<=5.0:
            saDict[ts.time].append(1)
        else:
            saDict[ts.time].append(0)

saList=[0]*len(resList)
for i in range(len(resList)):
    for k in saDict.keys():
        saList[i]+=saDict[k][i]

saList=[x/len(saDict.keys()) for x in saList]
saListno90=runningMean(saList,20)
resListno90=runningMean(resList,20)

print('done no 90\n')

u=MDAnalysis.Universe('NoGlycan/md.tpr','NoGlycan/md.xtc')
prot=u.select_atoms('protein')
sol=u.select_atoms('resname WN or resname CL or resname NA')

saDict={}
for ts in u.trajectory:
    saDict[ts.time]=[]
    for res in prot.residues[:-1]:
        minD=np.min(cdist(res.atoms.positions,sol.atoms.positions))
        if minD<=5.0:
            saDict[ts.time].append(1)
        else:
            saDict[ts.time].append(0)


saList=[0]*len(resList)
for i in range(len(resList)):
    for k in saDict.keys():
        saList[i]+=saDict[k][i]

saList=[x/len(saDict.keys()) for x in saList]
saListnoglycan=runningMean(saList,20)

plt.figure()
plt.plot(resListWT,saListWT,c='xkcd:orange',alpha=0.5,label='WT')
plt.plot(resListWT,saListno90,c='xkcd:red',alpha=0.5,label='-N90')
plt.plot(resListWT,saListnoglycan,c='xkcd:grey',alpha=0.5,label='deglycosylated')
plt.xlabel('residue')
plt.ylabel('proportion of time solvent-exposed')
plt.xticks([1,100,200,300,400,500,600,700],fontsize=8)
plt.legend()
plt.tight_layout()
plt.savefig('sasa.png',dpi=600)
plt.close()
