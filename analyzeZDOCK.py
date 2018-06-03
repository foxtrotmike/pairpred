# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:51:06 2013
Plot the comparison between PAIRpred and ZDOCK
@author: root
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import myPickle as mPickle
from DBD4 import parseCSVData
def getSortedR(cid,bdir,A,N=2000,M=2000):
    ifname=os.path.join(bdir,cid+'.zd3.0.2.cg.out.rmsds')
    R=np.zeros(N)
    with open(ifname, 'r') as f:
        for n in range(N):
            R[n]=f.readline().split()[1]
    sidx=np.argsort(A[cid][0][:M])
    RS=R+0.0
    
    RS[:M]=R[:M][sidx]
    return RS,R

bdir='../zdock_data/decoys_bm4_zd3.0.2_15deg/results/'
dfname='../ascores_2K2.scr.pkl'
d4=parseCSVData('..\Complete Data\DBD4_data.csv')
A=mPickle.load(dfname)
RS=[]
R=[]
dthr=2.5
N=2000
ncids=0
for i,cid in enumerate(A.keys()):
#    if cid not in d4 or d4[cid][1]=='RB':
#        continue
    rs,r=getSortedR(cid,bdir,A,N=N)
#    if np.any(r[:10]<dthr):
#        # import pdb
#        print cid, d4[cid][1]
#        
#        pdb.set_trace()
    R.append(np.cumsum(r<dthr)>0)
    RS.append(np.cumsum(rs<dthr)>0)
    ncids=ncids+1.0
plt.plot(range(1,N+1),np.sum(np.array(R),axis=0)/ncids,'k.-',label='ZDOCK')        
plt.plot(range(1,N+1),np.sum(np.array(RS),axis=0)/ncids,'ro-',label='ZDOCK with Resorting')    
plt.xlabel('Number of Predictions')    
plt.ylabel('Success Rate')
plt.axis([0,10,0,1])
plt.legend(loc=0)
plt.grid()
plt.show()