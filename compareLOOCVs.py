# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 17:16:13 2013

@author: root
"""
import glob
import numpy as np
from analyzeLOOCV import computeNTP, calcRFPP
dir1='./DBD3LOOCVSR/'
dir2='./DBD3LOOCV/'
fs=glob.glob(dir1+'*.pairpred.txt')
from BISEPutils import getFileParts
F1=np.zeros((len(fs),5))
F2=np.zeros((len(fs),5))
DNTP1=[]
DNTP2=[]
for i,f1 in enumerate(fs):
    
    fp=getFileParts(f1)
    f2=dir2+fp[1]+fp[2]
    (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx)=computeNTP(f1);
    DNTP1.append(dntp)
    print i,f1,auc, ttp, fpi,la,ra
    F1[i,:]=[auc ,ttp ,fpi ,la ,ra]
    (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx)=computeNTP(f2);
    DNTP2.append(dntp)
    print i,f2,auc, ttp, fpi,la,ra
    F2[i,:]=[auc ,ttp ,fpi ,la ,ra]
    
FPI1=F1[:,2]
FPI2=F2[:,2]
print "RFPP for",dir1
calcRFPP(FPI1,DNTP1)
print "RFPP for",dir2
calcRFPP(FPI2,DNTP2)

"""
RFPP for ./DBD3LOOCVSR/
[1.0, 4.0, 21.5, 60.0, 293.7000000000001]
[1.0, 2.75, 10.0, 39.25, 178.50000000000011]
[1.0, 2.0, 7.0, 27.0, 80.800000000000011]
[1.0, 1.0, 6.0, 17.25, 72.600000000000023]
[1.0, 1.0, 6.0, 16.25, 62.80000000000004]
RFPP for ./DBD3LOOCV/
[0.30000000000000071, 4.0, 20.5, 79.75, 275.50000000000017]
[1.0, 2.75, 11.0, 46.75, 172.69999999999999]
[1.0, 2.0, 7.0, 35.0, 128.50000000000006]
[1.0, 1.0, 5.0, 22.25, 104.80000000000001]
[1.0, 1.0, 3.0, 18.0, 71.400000000000006]

RFPP for ./DBD3LOOCVSR/
[1.0, 4.0, 20.5, 56.75, 275.10000000000025]
[1.0, 2.25, 10.0, 39.75, 157.50000000000028]
[1.0, 2.0, 7.0, 27.0, 78.400000000000034]
[1.0, 1.0, 6.0, 17.75, 67.800000000000068]
[1.0, 1.0, 5.5, 16.75, 54.400000000000119]
RFPP for ./DBD3LOOCV/
[0.0, 4.0, 20.0, 74.25, 242.50000000000045]
[1.0, 2.0, 10.0, 48.25, 172.10000000000002]
[1.0, 2.0, 7.0, 37.0, 119.50000000000013]
[1.0, 1.0, 4.5, 22.0, 102.40000000000003]
[1.0, 1.0, 3.0, 17.0, 70.200000000000017]

"""
