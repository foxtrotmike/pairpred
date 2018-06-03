# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 08:40:33 2013

@author: root
"""
from myPDB import *
from getExamplesDBD_breakup import *
from PyML.evaluators.roc import roc

E=getExamplesDBD.loader(os.path.join('../../DBD4CSPKL/PKL','ENS_15_35_50.lbl.pkl'))
pdbdir='../../DBD4CSPKL/PDB_all_'
pkldir='../../DBD4CSPKL/PKL'
F=list(set([getFileParts(g)[1].split('.')[0] for g in glob.glob(os.path.join(pkldir,'*.pdb.pkl'))]))
A={}
for fid in F:    
    print fid
    X=myPDB.loader(os.path.join(pkldir,fid+'.pdb.pkl'))
    C=np.max(X.pssm,axis=0)
    #C=X.rasa#np.sum(.psfm,axis=0)
    #C=JSON2ConsScore(ipdbfile, jfile)        
    fcids=[k for k in E.Pex.keys() if (fid in k)]
    fPi=[]
    for c in fcids:        
        fPi.extend([i[int(c[0]!=fid)] for i in E.Pex[c][0]])
    fPi=np.unique(np.array(fPi))
    if len(fPi):
        
        L=np.zeros(len(C))
        L[fPi]=1.0
        A[fid]=roc(list(C),list(L))[-1]
