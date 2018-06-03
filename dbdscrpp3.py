# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 23:32:18 2012

@author: fayyaz
Description: Performs cross validation on a kernel file (reduced data, not leave one complex out CV)
"""
from dbdKernels import dbdKernel, getKernelData
from myKernelData import *
from BISEPutils import *
from myPyMLassess import *
from PyML import *
from PyML.evaluators.assess import trainTest
from mpi4py import MPI
import numpy as np
import sys
from time import time
from PyML.evaluators import roc     
import random
import itertools
def getComplexFoldList(mK,nfolds=4,shuffle=False):    
    """
    If nfolds is a list of lists containing the complex ids for each fold then that is used
    otherwise it is created
    """
    dkey={}    
    for i,e in enumerate(mK.labels.patternID):
        try:
            v=dkey[e[0]]
        except:
            dkey[e[0]]=([],[])            
        dkey[e[0]][0].append(i)
        dkey[e[0]][1].append(e)
    
    
    if type(nfolds) is type([]):
        fcids=nfolds
        nfolds=len(fcids)
        #the fcids and dkey must have exactly the same elements
        vkeys=list(set(itertools.chain(*fcids)).intersection(dkey.keys()))
        dkey2={}
        for k in vkeys:
            dkey2[k]=dkey[k]
        dkey=dkey2
    else:
        cids=dkey.keys()
        if type(cids[0])!=type(''):
            cids0=list(set([c[0][:4] for c in cids]))
        else:
            cids0=cids
        if shuffle:
            random.shuffle(cids0)
        fcids0=list(chunks(cids0,int(np.ceil(len(cids0)/float(nfolds)))))
        
        if type(cids[0])!=type(''):
            fcids=[]
            for fx in fcids0:
                fcids.append([c for c in cids if c[0][:4] in fx])
        else:
            fcids=fcids0
    nfolds=len(fcids)
    print "Number of actual folds",nfolds
    print fcids
    
    F=[[] for f in range(nfolds)]
    for f in range(nfolds):
        
        for cid in fcids[f]:
            if cid in dkey:
                F[f].extend(dkey[cid][0])
            else:
                print "Warning: Complex",cid,"Not found in the kernel. Will not be used in the analysis"
    return F,dkey

def getAUC(s):
    if type(s)==type(''):
        (r,dkey)=cPickle.load(open(s, "rb" ) )
    else:
        (r,dkey)=s
 
    patid=combineList(r.getPatternID())
    vkey=dict(zip(patid,range(len(patid))))
    decfn=combineList(r.getDecisionFunction())
    lblid=combineList(r.getGivenLabels())
    cids=dkey.keys()
    D=[[] for i in cids]
    L=[[] for i in cids]
    A=[[] for i in cids]
    try:
        R=getRMSDDict('shandar_rmsd.txt')
    except:
        R=None
    Rx=[[] for i in cids]
    for i,cid in enumerate(cids):
        cidx=dkey[cid]        
        if type(cidx) is tuple: #backward compatability to old results objects 
            cidx=cidx[0]
        for e in cidx:
            try:
                n=vkey[e]
            except KeyError:
                pdb.set_trace()
            D[i].append(decfn[n])
            L[i].append(lblid[n])
        (_,_,a)=roc.roc(D[i],L[i])
        A[i]=a
        if R is not None:
            Rx[i]=R[cid]        
    (fp,tp,auc)=roc.roc_VA(zip(D,L))
    return (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))
def computeCvals(S,E):
    # first assign C based on the distance between residues in an example 
    lowC=0.01
    pdthr=E.pdthr
    ww=(E.mdthr-E.pdthr)/4.0 # 0.5
    gdf=lambda x: 1-(1-lowC)*np.exp(-((x-pdthr)**2)/(2*(ww)**2)) # the slack rescaling function
    Cintra=[+1.0 for _ in S]
    for i,(cid,e,_) in enumerate(S):
        if cid in E.dict_gdf and e in E.dict_gdf[cid]:
            Cintra[i]=gdf(E.dict_gdf[cid][e])       
    # Assign C based on distance of an example from its nearest example of opposite sign
    dpdlist=list(itertools.chain(*[E.dpd[c].values() for c in E.dpd.keys()]))
    dndlist=list(itertools.chain(*[E.dnd[c].values() for c in E.dnd.keys()]))
    ppmin=np.min(dpdlist)
    npmin=np.min(dndlist)
    pp25=np.percentile(dpdlist,10)
    np25=np.percentile(dndlist,10)
    ldf=lambda x,xmin,xmax: lowC+((x-xmin)*(1-lowC)/(xmax-xmin))
    Cinter=[+1.0 for _ in S]
    for i,(cid,e,_) in enumerate(S):
        if cid in E.dpd and e in E.dpd[cid]:
            f=ldf(E.dpd[cid][e],ppmin,pp25)
            #print f
            if f<1.0:
                Cinter[i]=f
        elif cid in E.dnd and e in E.dnd[cid]:
            f=ldf(E.dnd[cid][e],npmin,np25)
            if f<1.0:
                Cinter[i]=f
    for i in range(len(Cinter)): #take the min
        if Cintra[i]<Cinter[i]:
            Cinter[i]=Cintra[i]
    return Cinter
            
def getRMSDDict(fname):
    D={}
    for l in open(fname,'r'):
        x=l.split()
        D[x[0][:4]]=float(x[1])
    return D
    
if __name__=="__main__":      
    pwKfname='dbdKernel.dbk.pkl'
    print pwKfname
    s=SVM()
    s.C=10   
    print "Loading Kernel from file..."
    t1=time()
    mK=getKernelData(*dbdKernel.loader(pwKfname)) #NOTE THAT OBJECT IS BEING LOADED FROM FILE
    y=mK.labels.L    
    s.Cmode='classProb'
    if s.Cmode=='fromData':
        s.Cmode='classProb';Clist=s.getClist(mK);s.Cmode='fromData'
        (S,_,E)=dbdKernel.loadNonKernelData(pwKfname)  
        Ck=computeCvals(S,E)        
        C=[Clist[i]*Ck[i] for i in range(len(Ck))]
        mK.C=C
        mK.registerAttribute('C',mK.C)
    #pdb.set_trace()
    t2=time()
    print "Pairwise kernel loaded in",t2-t1,"s"   
    
    ttP,dkey=getComplexFoldList(mK)
    
    mK.attachLabels(Labels(y))    #pyml did not work with the string labels
    print "Starting CV ... "
    #r=mycvFromFolds(s,mK,testingPatterns=ttP)#,fcnTrainTest=pegasostrainTest) 
    #"""
    #from pegasos import *
    #s=Pegasos()
    C=10.0
    s.Lambda=1.0/((len(mK)/len(ttP))*C)
    s.nCycles=30000#np.int(100/s.Lambda)
    print s
    r=mycvFromFolds(s,mK,testingPatterns=ttP) 
    #"""
    if r is not None:
        output = open('result.res.pkl', 'wb')
        cPickle.dump((r,dkey), output)
        output.close()
        (auc,_,_)=getAUC((r,dkey))
        print "CV Complete. Total Time taken (s):",time()-t1
        print "Complex-wise average AUC",auc
        print "Overall AUC = ", r.getROC()