# -*- coding: utf-8 -*-
"""
Created on Mon May  6 20:40:10 2013
DI correlation analysis: Analyze the correlation between inter-residue distance and MI and DI
@author: root
"""
import os
import numpy as np
# import pdb
#from BISEPutils import getFileParts
from hhblits import runHHblits
from calcEC import runEC, parseECFile
from myPDB import *
def toArray(L,d):
    Dm=np.zeros((len(L.Coords),len(L.Coords)))
    Dm.fill(np.nan)
    for k in d:
        Dm[k[0],k[1]]=d[k]
    return Dm
def getDist(L):
    D=np.zeros((len(L.Coords),len(L.Coords)))
    D.fill(np.nan)
    for i in range(len(L.Coords)):
        for j in range(i,len(L.Coords)):
            d=spatial.distance.cdist(L.Coords[i], L.Coords[j]).min()
            D[i,j]=d#L.rASA[i]*L.rASA[j]
    return D
def getECForProtein(pdbpklpath,fid,odir='./'):
    """
    given a complex 'cid', it computes the MSA and MI and DI
    """
    # Write the ligand
    L=myPDB.loader(os.path.join(pdbpklpath,fid+'.pdb.pkl'))    
    lfasta=os.path.join(odir,fid+'.seq')    
    lafasta=lfasta+'.fasta'    
    ecfile=os.path.join(odir,fid+'.ec')    
    if os.path.exists(ecfile) and os.path.isfile(ecfile):#If the EC file already exists
        print "Using existing EC file", fid, ecfile
        (M,D)=parseECFile(ecfile) 
    else:
        if not (os.path.exists(lafasta) and os.path.isfile(lafasta)):#if the file lafasta does not exist
            print "Calculating MSA for", lafasta
            L.save2FASTA(lfasta,saveHeader=True,hdr=fid[:4]+'/'+fid[4:])  
            runHHblits(lfasta,niter=2,cpu=1)
        else:
            print "Using existing MSA for", lafasta
        print "Calculating EC for", fid
        (M,D)=runEC([lafasta],fid[:4],ofile=ecfile)
    Mi=toArray(L,M).flatten()
    Di=toArray(L,D).flatten()
    Dl=getDist(L).flatten()
    nanidx=~(np.isnan(Mi)+np.isnan(Di)+np.isnan(Dl))
    Mi=Mi[nanidx]
    Di=Di[nanidx]
    Dl=Dl[nanidx]
    return Mi,Di,Dl
    #get the seq for L and R and use it to get the 
def getAUCs(Mx,Dx,D,dthr=6.0):
    L=list(2*(D<dthr)-1)
    (_,_,aa_di)=roc.roc(list(Dx),L)        
    (_,_,aa_mi)=roc.roc(list(Mx),L)
    return aa_di,aa_mi
if __name__=="__main__":

    pdbpklpath='./DBD4N/PDBPKL4'#'/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'
    evalROC=True
    if evalROC:
        try:
            from PyML.evaluators import roc        
        except ImportError:
            print "PyML not found, will not evaluate ROC."
            evalROC=False 
    """       
    fid='1SBB_l_u'
    Mx,Dx,D=getECForProtein(pdbpklpath,fid,odir='./')
    if evalROC:
        aa_di,aa_mi=getAUCs(Mx,Dx,D,dthr=6.0)
        print 'AUC_DI=',aa_di
        print 'AUC_MI=',aa_mi
#    print np.corrcoef(Mx,D)
#    print np.corrcoef(Dx,D)
     """
    
        
    incids=['2I25', '2MTA', '2NZ8', '2O3B', '2O8V', '2OOR', '2OZA', '2SNI', '3D5S', '3SGQ', 'BOYV', '1QFW', '1R0R', '1SBB', '1T6B', '1US7', '1WDW', '1WEJ', '1Z0K', '2BTF', '2CFH', '2H7V', '1IQD', '1JMO', '1JWH', '1K4C', '1KKL', '1KXP', '1M10', '1MLC', '1MQ8', '1N8O', '1NSN', '1OYV', '1FSK', '1H9D', '1HCF', '1HE1', '1HE8', '1HIA', '1I2M', '1I4D', '1IJK', '1ACB', '1AK4', '1AY7', '1BJ1', '1BKD', '1DQJ', '1E96', '1EER', '1F34', '1F6M', '1FC2', '1FCC']
    odir='.\ecstats_single'
    AUC_DI=[]
    AUC_MI=[]
    for f in incids:
        print "processing",f
        fid=f+'_l_u'
        Mx,Dx,D=getECForProtein(pdbpklpath,fid,odir=odir)
        aa_di,aa_mi=getAUCs(Mx,Dx,D)
        AUC_DI.append(aa_di)
        AUC_MI.append(aa_mi)
        fid=f+'_r_u'
        Mx,Dx,D=getECForProtein(pdbpklpath,fid,odir=odir)
        aa_di,aa_mi=getAUCs(Mx,Dx,D)
        AUC_DI.append(aa_di)
        AUC_MI.append(aa_mi)        
    
    print np.mean(AUC_DI), np.mean(AUC_MI)
    