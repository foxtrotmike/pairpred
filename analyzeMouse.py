# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 19:46:32 2013
Code for finding out the differences in prediction scores for the mouse ISG15 in case study of PAIRpred paper
@author: root
"""
import numpy as np
from analyzePredFile import *
import matplotlib.pyplot as plt
mppfile='./testcase/DBD4_IGS15_NS1B_PairPred/Mutagenesis/M6B1_D4.InterPRed.txt'
hppfile='./testcase/DBD4_IGS15_NS1B_PairPred/HNS1_d4PN.InterPRed.txt'#HNS1.InterPRed200.txt'
(_,hMv,_,hlseq,_,_,_)=readFile(hppfile)
(_,mMv,_,mlseq,_,_,_)=readFile(mppfile)

"""
nz=np.zeros((1,mMv.shape[1]))
nz.fill(np.nan)
mMv=np.vstack((mMv[:48,:],nz,mMv[47:51,:],nz,mMv[51:,:]))[:hMv.shape[0],:]
mMv[:9,:]=np.nan;mMv[-6:,:]=np.nan;mMv[:,:7]=np.nan;mMv[:,-6:]=np.nan;
hMv[:9,:]=np.nan;hMv[-6:,:]=np.nan;hMv[:,:7]=np.nan;hMv[:,-6:]=np.nan;
hv=np.nanmax(hMv,axis=1)
mv=np.nanmax(mMv,axis=1)
mv=(mv-np.mean(mv[~np.isnan(mv)]))/np.std(mv[~np.isnan(mv)])
hv=(hv-np.mean(hv[~np.isnan(hv)]))/np.std(hv[~np.isnan(hv)])
print np.nonzero((hv>1.3) & ((mv-hv)<-0.8))
plt.plot(hv);plt.plot(mv);plt.show()
"""
nz=np.zeros((1,mMv.shape[1]))
nz.fill(np.nan)
mMv=np.vstack((mMv[:48,:],nz,mMv[47:51,:],nz,mMv[51:,:]))[:hMv.shape[0],:]
#mMv[:9,:]=np.nan;mMv[-5:,:]=np.nan;mMv[:,:6]=np.nan;mMv[:,-5:]=np.nan;
#hMv[:9,:]=np.nan;hMv[-5:,:]=np.nan;hMv[:,:6]=np.nan;hMv[:,-5:]=np.nan;
mMv=(mMv-np.mean(mMv[~np.isnan(mMv)]))/np.std(mMv[~np.isnan(mMv)])
hMv=(hMv-np.mean(hMv[~np.isnan(hMv)]))/np.std(hMv[~np.isnan(hMv)])
hMvf=hMv.flatten()
hMvf=hMvf[~np.isnan(hMvf)]
thr=-10.0
dthr=1.0
R=mMv-hMv
Rf=R.flatten()
Rf=Rf[~np.isnan(Rf)]
R[hMv<=thr]=np.nan
hMv[hMv<=thr]=np.nan
R[R>-dthr]=np.nan
#R[mMv>thr-dthr]=np.nan
print np.unique(np.nonzero(~np.isnan(R))[0])
print np.nonzero(~np.isnan(np.nanmax(R,axis=1)))
dd=dict(zip(*(np.nonzero(~np.isnan(np.nanmax(R,axis=1)))[0],np.nanmax(R,axis=1)[np.nonzero(~np.isnan(np.nanmax(R,axis=1)))[0]])))
print dd
#plt.figure();plt.hist(hMvf,cumulative=True,normed=True,bins=200)


#plt.figure();plt.hist(Rf,cumulative=True,normed=True,bins=50)
p=(np.sum(hMv>thr)/float(hMv.size))*100
print p
"""
cmap =plt.cm.jet;# plt.get_cmap()
cmap.set_bad(color = 'w', alpha = 1.)

plt.figure();plt.imshow(R,cmap=cmap);plt.colorbar();
plt.title('Considering top %1.1f predictions with decrease in prediction score' %p)
plt.show()
"""