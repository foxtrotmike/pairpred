# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 16:39:23 2013

@author: root
"""
from myPDB import *
import random
#Build constraint vector
#def buildConstraints(L,R,**args):


def buildConstraints(L,R,Nc=1000):
    DL=getDistMat(L.Coords); DR=getDistMat(R.Coords); 
    DLT=np.median(DL);DRT=np.median(DR);
    DL=DL/DLT;DR=DR/DRT
    DT=2.0
    THR_rasa=0;
    LD=(DL>DT)*(L.rASA>THR_rasa)
    RD=DR>DT*(R.rASA>THR_rasa)
    vLi=np.nonzero(np.any(LD,axis=0))[0]
    vRi=np.nonzero(np.any(RD,axis=0))[0]
    E=[]
    n=0
    ntries=100    
    while (n<Nc):
        l1=random.sample(vLi,1)[0]
        r1=random.sample(vRi,1)[0]
        l2=random.sample(np.nonzero(LD[l1,:])[0],1)[0]
        r2=random.sample(np.nonzero(RD[r1,:])[0],1)[0]    
        e=((l1,r1),(l2,r2))
        e_=((l1,r2),(l2,r1))
        if e not in E and tuple(reversed(e)) not in E and e_ not in E and tuple(reversed(e_)) not in E:
            E.append(e)            
            n+=1
            ntries=100
        ntries=ntries-1
        if ntries<=0:
            print "Warning: Required number of constraints could not be made. Made: ", len(E)
            break;
    return E
# Find the example number in test data of (l1,r1) and (l2,r2)

if __name__=="__main__":
    cid='1A2K'
    bdir='../DBD4N/PDBPKL4/'
    L=myPDB.loader(bdir+cid+'_l_u.pdb.pkl')
    R=myPDB.loader(bdir+cid+'_r_u.pdb.pkl')


    