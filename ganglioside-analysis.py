#!/usr/bin/python

import numpy as np
import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist
#25,1-4; 3 is SIA
def runningMean(x,N):
    cumsum=np.cumsum(np.insert(x,0,0))
    return (cumsum[N:]-cumsum[:-N])/float(N)

sns.set_style('white')

u=MDAnalysis.Universe('md.tpr','md.xtc')
phospholipid=u.select_atoms('resname POPC')
ganglioside=u.select_atoms('resname HXS or resname SIA or resname CER')
L=LeafletFinder(u,'name NC3')
leaflet0=L.groups(0)
leaflet1=L.groups(1)

heightDict={}
datDict={
        (None,4):[],
        (4,8):[],
        (8,14):[]
        }
mapDict={
        (None,4):'HXS1',
        (4,8):'HXS2',
        (8,14):'SIA'
        }
for ts in u.trajectory:
    heightDict[ts.time]=[]

splitgangliosides=[]
j=None
for i in np.arange(25,25*20,25):
    splitgangliosides.append(ganglioside.atoms[j:i])
    j=i

timeList=[]
hgContacts=[]
amContacts=[]
for ts in u.trajectory:
    timeList.append(ts.time)
    seen=[]
    count1=0
    count2=0
    for ag1 in splitgangliosides:
        for ag2 in splitgangliosides:
            if ag1==ag2:
                next
            elif (ag1,ag2) in seen or (ag2,ag1) in seen:
                next
            else:
                seen.append((ag1,ag2))
                if np.min(cdist(ag1.atoms[:15].positions,ag2.atoms[:15].positions))<=6.0:
                    count1+=1
                if np.min(cdist(ag1.atoms[16:18].positions,ag2.atoms[16:18].positions))<=6.0:
                    count2+=1
    hgContacts.append(count1)
    amContacts.append(count2)

plt.figure()
sns.set_style('white')
plt.plot(timeList,hgContacts,c='xkcd:grey',alpha=0.5,lw=1.5)
plt.plot(timeList,amContacts,c='xkcd:grey',alpha=0.5,lw=1.5)
plt.plot(runningMean(timeList,50),runningMean(hgContacts,50),c='xkcd:red',label='glycans')
plt.plot(runningMean(timeList,50),runningMean(amContacts,50),c='xkcd:orange',label='AM1')
plt.xlabel('time (ps)')
plt.ylabel('number of contacts')
plt.legend()
plt.tight_layout()
plt.savefig('gm3-contacts.png',dpi=600)
plt.close()

for ts in u.trajectory:
    area=u.dimensions[0]*u.dimensions[1]
    heightList=[]
    for gm1 in splitgangliosides:
        hgpos=gm1.atoms[16:18].center_of_geometry()[2]
        if np.linalg.norm(hgpos-leaflet0.atoms.center_of_geometry()[2])<np.linalg.norm(hgpos-leaflet1.atoms.center_of_geometry()[2]):
            leaflet=leaflet0
        else:
            leaflet=leaflet1
        tmp=[]
        for pair in [(None,4),(4,8),(8,14)]:
            ag=gm1[pair[0]:pair[1]]
            height=np.linalg.norm(ag.atoms.center_of_geometry()[2]-leaflet.atoms.center_of_geometry()[2])
            datDict[pair].append(np.linalg.norm(ag.atoms.center_of_geometry()[2]-leaflet1.atoms.center_of_geometry()[2]))
            tmp.append(height)
        heightList.append(np.max(tmp))
    heightDict[ts.time]=heightList

datList=[]
for k in heightDict.keys():
    datList.append(np.mean(heightDict[k]))
plt.figure()
sns.distplot(datList)
plt.xlabel(r'maximum height ($\AA$)')
plt.ylabel('frequency')
plt.tight_layout()
plt.savefig('gm3-heights.png',dpi=600)
plt.close()

plt.figure()
datList=[]
stdList=[]
xList=[]
nameList=[]
for pair in datDict.keys():
    xList.append(pair)
    datList.append(np.mean(datDict[pair]))
    stdList.append(np.std(datDict[pair]))
plt.scatter(range(len(datList)),datList,c='xkcd:orange')
plt.errorbar(range(len(datList)),datList,yerr=stdList,c='xkcd:grey',ls='None',alpha=0.5)
for x in xList:
    nameList.append(mapDict[x])
plt.xticks(range(len(xList)),nameList)
plt.xlabel('group')
plt.ylabel(r'mean height ($\AA$)')
plt.tight_layout()
plt.savefig('gm3-groupheights.png',dpi=600)
plt.close()

L=LeafletFinder(u,'name AM1 or name PO4')

timeList=[]
aplPOPC=[]
for ts in u.trajectory[:10]:
    area=u.dimensions[0]*u.dimensions[1]
    aplPOPC.append(area/len(L.groups(1)))
aplPOPC=np.mean(aplPOPC)
print(aplPOPC)
aplDPG3=[]
aplList=[]
for ts in u.trajectory:
    timeList.append(ts.time)
    tmp=[]
    area=u.dimensions[0]*u.dimensions[1]
    phosCount=0
    for a in L.groups(0).atoms.names:
        if a=='PO4':
            phosCount+=1
    aplDPG3.append((area-(phosCount*aplPOPC))/(len(L.groups(0).atoms)-phosCount))
    tmp.append((area-(phosCount*aplPOPC))/(len(L.groups(0).atoms)-phosCount))
    aplList.append(np.mean(tmp))

plt.figure()
plt.plot(timeList,aplList,c='xkcd:grey',alpha=0.5,lw=1.5)
plt.plot(runningMean(timeList,50),runningMean(aplList,50),c='xkcd:orange',label='%s %s'%(np.mean(aplDPG3),np.std(aplDPG3)))
plt.xlabel('time (ps)')
plt.ylabel(r'area per lipid ($\AA^{2}$)')
plt.legend()
plt.tight_layout()
plt.savefig('gm3-apl.png',dpi=600)
plt.close()
