# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:47:23 2013

@author: root
"""

from dbdscrpp3 import *
import itertools
import matplotlib.pyplot as plt

# import pdb
if __name__=="__main__":
    fname='../Results/result_tppk.res.pkl'
    pdbpklpath='../PDBPKLP3'
    (auc,(fp,tp),(A,Rx,Dx,Lx,cids,r,dkey))=getAUC(fname)    
    Mv=np.zeros(len(cids))
    Mv.fill(np.nan)
    Mvtpr=np.zeros(len(cids))
    Mvtpr.fill(np.nan)
    for (i,cid) in enumerate(cids):
        d=Dx[i]
        l=np.array(Lx[i])
        (c,b,p)=plt.hist(d,normed=True,cumulative=True,bins=50)
        b=b[::-1]
        b=b[1:]
        c=1-c[::-1]
        for ib,bv in enumerate(b):
            if np.any(l[d>bv]==-1):
                break
        Mvtpr[i]=100*np.sum(l[d>bv]==+1)/float(np.sum(l==+1))
        Mv[i]=100*c[ib]
        #pdb.set_trace()
    """
    D=np.array(list(itertools.chain(*Dx)))
    L=np.array(list(itertools.chain(*Lx)))
    
    #PLOTING THE TPR and FPR versus the threshold
    plt.figure(0);(c,b,p)=plt.hist(D,normed=True,cumulative=True,bins=100)
    #mD=np.nanmin(D)
    #MD=np.nanmax(D)
    b=b[:-1]
    #b=np.linspace(mD,MD,100)
    fpr=np.zeros(c.shape)
    fpr.fill(np.nan)
    tpr=np.zeros(c.shape)
    tpr.fill(np.nan)
    ntp=float(np.sum(L==+1))
    ntn=float(np.sum(L==-1))
    for (i,bv) in enumerate(b):
        pidx=D>bv
        nidx=D<bv
        
        fpr[i]=np.sum(L[pidx]==-1)/ntn
        tpr[i]=np.sum(L[pidx]==+1)/ntp
    plt.figure(1);plt.plot(1-c,fpr,'r',1-c,tpr,'b');
    
    plt.show()
    """