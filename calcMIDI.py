# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:17:59 2013

@author: root
"""
from PyML.evaluators import roc       
from myPDB import *
import os
from hhblits import runHHblits
from calcEC import runEC, parseECFile
import numpy as np
# import pdb

def getDist(LCoords,RCoords=None):
    sym=False
    if RCoords is None:
        RCoords=LCoords
        sym=True
    D=np.zeros((len(LCoords),len(RCoords)))
    D.fill(np.nan)
    stidx=0
    for i in range(len(LCoords)):
        if sym:
            stidx=i       
        for j in range(stidx,len(RCoords)):
            d=spatial.distance.cdist(LCoords[i], RCoords[j]).min()
            D[i,j]=d#L.rASA[i]*L.rASA[j]
    return D
def toArray(d,ll,lr=None):
    if lr is None:
        lr=ll
    Dm=np.zeros((ll,lr))
    Dm.fill(np.nan)
    for k in d:
        Dm[k[0],k[1]]=d[k]
    return Dm
def getAUCs(mx,dx,dd,dthr=6.0):
    lx=2*(dd<dthr)-1
    mxf=mx.flatten()
    dxf=dx.flatten()
    lxf=lx.flatten()
    nanidx=~(np.isnan(mxf)+np.isnan(lxf)+np.isnan(dxf))
    lxf=list(lxf[nanidx])
    mxf=list(mxf[nanidx])
    dxf=list(dxf[nanidx])
    (_,_,aa_mi)=roc.roc(mxf,lxf)
    (_,_,aa_di)=roc.roc(dxf,lxf)
    return aa_mi,aa_di
def getECFromFile(fname,ll,lr):
    Mx={}
    Dx={}
    (M,D)=parseECFile(fname)
    for l in range(ll):        
        for r in range(lr):
            Mx[(l,r)]=M[(l,r+ll)]
            Dx[(l,r)]=D[(l,r+ll)]
    return Mx,Dx
    
def getECForComplex(pdbpklpath,cid,odir='./'):
    """
    given a complex 'cid', it computes the MSA and MI and DI
    """
    # Write the ligand
    L=myPDB.loader(os.path.join(pdbpklpath,cid+'_l_u.pdb.pkl'))
    R=myPDB.loader(os.path.join(pdbpklpath,cid+'_r_u.pdb.pkl'))
    #pdb.set_trace()
    lfasta=os.path.join(odir,cid+'_l_u.seq')
    rfasta=os.path.join(odir,cid+'_r_u.seq')
    lafasta=lfasta+'.fasta'
    rafasta=rfasta+'.fasta'
    ecfile=os.path.join(odir,cid+'.ec')
    if os.path.exists(ecfile) and os.path.isfile(ecfile):#If the EC file already exists
        print "Using existing EC file", cid, ecfile
        (M,D)=parseECFile(ecfile) 
    else:
        if not (os.path.exists(lafasta) and os.path.isfile(lafasta)):#if the file lafasta does not exist
            print "Calculating MSA for", lafasta
            L.save2FASTA(lfasta,saveHeader=True,hdr=cid+'/l')  
            runHHblits(lfasta,niter=2,cpu=1)
        else:
            print "Using existing MSA for", lafasta
        if not (os.path.exists(rafasta) and os.path.isfile(rafasta)):#if the file rafasta does not exist
            print "Calculating MSA for", rafasta
            R.save2FASTA(rfasta,saveHeader=True,hdr=cid+'/r')    
            runHHblits(rfasta,niter=2,cpu=1)
        else:
            print "Using existing MSA for", rafasta
        print "Calculating EC for", cid
        (M,D)=runEC([lafasta,rafasta],cid,ofile=ecfile)
    ll=len(L.seq)
    lr=len(R.seq)
    #pdb.set_trace()
    Mx={}
    Dx={}
    Ml={}
    Dl={}
    Mr={}
    Dr={}
    for l1 in range(ll):    
        for l2 in range(l1+1,ll):
            Ml[(L.S2Ri[l1],L.S2Ri[l2])]=M[(l1,l2)]
            Dl[(L.S2Ri[l1],L.S2Ri[l2])]=D[(l1,l2)]
    for l1 in range(lr):    
        for l2 in range(l1+1,lr):
            Mr[(R.S2Ri[l1],R.S2Ri[l2])]=M[(l1+ll,l2+ll)]
            Dr[(R.S2Ri[l1],R.S2Ri[l2])]=D[(l1+ll,l2+ll)]    
    for l in range(ll):        
        for r in range(lr):
            Mx[(L.S2Ri[l],R.S2Ri[r])]=M[(l,r+ll)]
            Dx[(L.S2Ri[l],R.S2Ri[r])]=D[(l,r+ll)]
    #pdb.set_trace()
    return (Mx,Dx,getDist(L.Coords,R.Coords)),(toArray(Ml,len(L.Coords)),toArray(Dl,len(L.Coords)),getDist(L.Coords)),(toArray(Mr,len(R.Coords)),toArray(Dr,len(R.Coords)),getDist(R.Coords))
    #get the seq for L and R and use it to get the 

if __name__=="__main__":
    cid='1SBB'
    pdbpklpath='./DBD4N/PDBPKL4'#'/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'
    Mx,Dx=getECForComplex(pdbpklpath,cid,odir='./')