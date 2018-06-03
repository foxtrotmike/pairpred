# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 20:41:49 2013
Code to get the AUC score and the RFPP for a CAPRI target (or any other complex which has a single complex file rather than having separate bound ligand and receptor PDB files) from its pairpred prediction file.
@author: root
"""
from jointAnalysis import *
from PyML.evaluators import roc
from postProcess import postProcessAvg
#from getExamplesDBD import getPosex
from symmetryProcessing import *
#bdir='../CAPRI/'
bdir='../../g2mers/'
cid='1MLC'
#(_,_,P,_,_)=getPosex(bdir,cid)    #get positive examples    
P=getPosexFromPDB(bdir,cid,dthr=6.0)        # Handles symmetry in the complex
ppfile=bdir+cid+'.pairpred.txt'

(auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(ppfile,usePDBidx=False)
#auc0,Mvc0,Mv,Mlc,lseq,rseq,lrV0,lrV,rrV0,rrV=postProcessAvg(cid,bdir,bdir)
#
#Mv[:10,:]=np.nan
#Mv[-10:,:]=np.nan
#Mv[:,:10]=np.nan
#Mv[:,-10:]=np.nan

Mvtbl=np.zeros(Mv.shape)
for (i,j) in P:
    Mvtbl[i,j]=1.0
Mvr=Mv.ravel()
Mvtblr=Mvtbl.ravel()
nidx=(~np.isnan(Mvr))
Mvr=Mvr[nidx]
Mvtblr=Mvtblr[nidx]
(fpv,tpv,aucv)=roc.roc(list(Mvr),list(Mvtblr))
print cid,"AUC =",aucv, "RFPP =",np.argmax(Mvtblr[np.argsort(-Mvr)]==1)
(fpl,tpl,auc)=roc.roc(list(np.nanmax(Mv,axis=0)),list(np.nanmax(Mvtbl,axis=0))); print auc
(fpl,tpl,auc)=roc.roc(list(np.nanmax(Mv,axis=1)),list(np.nanmax(Mvtbl,axis=1))); print auc