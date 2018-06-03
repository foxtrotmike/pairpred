# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 00:30:57 2012
@author: Fayyaz Minhas
Difference from dbdKernels2.py: This file does not save the kernel, instead it 
directly produces the results
This module implementes different kernels for the dbd data. It allows constr-
uction of pairwise kernel given the just the pdbpkl files and the getExamplesDBD
object corresponding to them.

If you run the module then it will ecaluate the dbdKernel and save it
    mpirun -np n python dbdKernels2.py
When ussing the saved dbdKernel file 
    from dbdKernels import pwKernel,dbdKernel,complexKernel
"""
from BISEPutils import *
from myPDB import *
from getExamplesDBD import *
from protKernel import *
from myKernelData import *
import cPickle
from mpi4py import MPI
from copy import deepcopy
from time import time
from dbdscrpp3 import *
def computecKdict(KK,pdbpklpath,E,scid=None,comm=None,myid=0,nprocs=1):    
    """
    Computes the dictionary object cKdict
    """
    aA=False 
    if len(KK)==0:
        assert scid is not None
        KK=[scid]
        aA=True
    Nc=len(KK)     
    bA=False
    if scid is None:
        clist=[(a,b) for a in range(Nc) for b in range(a,Nc)]
    else:
        clist=[(a,scid) for a in range(Nc)]
        bA=True
    
    csize=int(np.ceil(len(clist)/float(nprocs)))
    gclist=list(chunks(clist,csize))
    myclist=gclist[myid]
    mycK={}
    for (cpa,cpb) in myclist:
        if bA:
            kcpb=cpb
        else:
            kcpb=KK[cpb]
        cpkey=(KK[cpa],kcpb)
        print "Processing complexes : ",cpkey
        mycK[cpkey]=complexKernel(KK[cpa],kcpb,pdbpklpath,E,aAll=aA,bAll=bA)
    
    cKdict=None  
    if(myid!=0):
        comm.send(mycK, dest=0)
        
    if(myid==0):
        gcK=[mycK]
        for p in range(1,nprocs):
            gcK.append(comm.recv(source=p))            
        #self.E=E
        #self.pdbpklpath=pdbpklpath
        cKdict=mergeDicts(gcK)
        
    return cKdict
    
def dbdK2pwK(cKdict,Sr,Sc=None,justdiag=False):
    sym=False
    if Sc is None:
        Sc=Sr
        sym=True
    if justdiag:
        assert sym
        Krc=np.zeros(len(Sr),np.single)    
    else:
        Krc=np.zeros((len(Sr),len(Sc)))#Krc=np.zeros((len(Sr),len(Sc)),np.single)
    st=0
    for aidx in xrange(len(Sr)):     
        #print aidx
        if justdiag:
            sa=Sr[aidx]
            tK=cKdict[(sa[0],sa[0])]
            Krc[aidx]=tK2kv(tK,sa,sa)
            continue
        if sym:
            st=aidx
        for bidx in xrange(st,len(Sc)):
            sa=Sr[aidx]              
            sb=Sc[bidx]
            #pdb.set_trace()
            try:
                tK=cKdict[(sa[0],sb[0])]
            except Exception as exc:
                if sym:
                    sa,sb=sb,sa
                    tK=cKdict[(sa[0],sb[0])]
                else:
                    raise exc            
            kv=tK2kv(tK,sa,sb)            
            Krc[aidx,bidx]=kv
            if sym:
                Krc[bidx,aidx]=kv
    return Krc
        
class dbdKernel:
    """
    Given a getExamplesDBD object, it evaluates the kernel values over all
    examples from all complexes in it
    E: the getExamplesDBD object
    pdbpklpath: path where the pdbpkl files are kept
    ofname: where the output file is stored (default: dbdKernel.dbK.pkl)
    cKdict: is a dictionary object which when given the (cidA,cidB) key will 
        return the complex kernel between them
        Since the kernel is evaluated between every two complexes so only the
        upper triangular matrix (along with the diagnoal) is stored
        so if you cannot find (cidA,cidB) try (cidB,cidA)
    cids: complex ids for which the kernel has been evaluated (obtained from E)
    Note: This object support MPI style parallelization and all you have to 
        do is to call the constructor in a program that has been run using 
        mpiexec or mpirun. Only processor 0 will contain the 'valid' (complete)
        object which is saved in ofname. Each process gets equal load.
    """
    def __init__(self,pdbpklpath,E0,ofname=None):
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        nprocs = comm.Get_size()
        
        E=getValidE(E0,pdbpklpath) 
        self.E=E
        KK=E.Pex.keys()
        if myid==0:
            txx0=time()
            ignored=list(set(E0.Pex.keys()).difference(KK))
            if(len(ignored)):
                print 'Ignoring the following complexes because at least one of the \
                four required pdbpkl files could not be found in', pdbpklpath,':'\
                ,ignored
            print "Evaluating kernel over: ", KK
            print "Using",nprocs,"processe(s)."
        #pdb.set_trace()
        cKdict=computecKdict(KK,pdbpklpath,E,comm=comm,myid=myid,nprocs=nprocs)
        if myid==0:            
            self.__dbdK2pwK__(KK,E,cKdict)
            """	
            print "TIME TAKEN FOR KERNEL EVALUATION (s): ",time()-txx0
            ttP,dkey=getComplexFoldList(self.mK)
            self.mK.attachLabels(Labels(self.mK.labels.L)) 
            t1=time()
            print "Starting CV ... "
            s=SVM()
            s.C=10 
            r=mycvFromFolds(s,self.mK,testingPatterns=ttP)
            if r is not None:
                output = open('result.res.pkl', 'wb')
                cPickle.dump((r,dkey), output)
                output.close()
                (auc,_,_)=getAUC((r,dkey))
                print "CV Complete. Total Time taken (s):",time()-t1
                print "Complex-wise average AUC",auc
                print "Overall AUC = ", r.getROC()
            """
            
            #f3=['1XQS','1A2K']
            #f4=['1ACB','1AHW','1XXX'] #'1A2K','1ACB','1AHW','1XQS
            f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
            f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
            s=SVM()
            s.C=10  
            t1=time()
            ttP,dkey=getComplexFoldList(self.mK,nfolds=[f3,f4])
            self.mK.attachLabels(Labels(self.mK.labels.L))    #pyml did not work with the string labels
            print "Starting CV ... "
            r=mycvFromFolds(s,self.mK,testingPatterns=ttP)
            if r is not None:
                output = open('result.res.pkl', 'wb')
                cPickle.dump((r,dkey), output)
                output.close()
                (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))=getAUC((r,dkey))        
                print "CV Complete. Total Time taken (s):",time()-t1
                print "AUC for testing on DBD 3:",np.mean([A[i] for i,c in enumerate(cids) if c in f3])
                print "AUC for testing on DBD 4:",np.mean([A[i] for i,c in enumerate(cids) if c in f4])
                print "Complex-wise AUC = ",auc
                print "Overall AUC = ", r.getROC()
            
    def __dbdK2pwK__(self,cids,E,cKdict):        
        S=[(cid,e,+1) for cid in cids for e in E.Pex[cid][0]]
        N=[(cid,e,-1) for cid in cids for e in E.Nex[cid]]
        S.extend(N)
        K=dbdK2pwK(cKdict,S)    
        #HANDLE NANS HERE: If a diagonal is a nan, remove it
        nanidx=~np.isnan(K.diagonal())        
        K=K[nanidx,:]
        K=K[:,nanidx]
        S=[s for (s,n) in zip(S,nanidx) if n]
        #Normalize
        dK=np.diag(K)
        K=K/np.sqrt(dK[:,np.newaxis]*dK[:,np.newaxis].T)        
        ##self.K=K
        ##self.S=S
        ##self.dK=dK
        self.mK=getKernelData(K,S)

    def save(self,ofname=None):
        """
        Save the object
        """
        if ofname is None:
            ofname='dbdKernel.dbK.pkl'
        output = open(ofname, 'wb')
        cPickle.dump((self.K,self.S,self.dK,self.E), output,-1)
        output.close()        
        print "File Saved: "+ofname
    @classmethod  
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """        
        return cPickle.load(open(pklfile, "rb" ) )
        
def tK2kv(tK,sa,sb):
    k_mlpk=(tK.K_ll[tK.lAidx[sa[1][0]],tK.lBidx[sb[1][0]]]+\
    tK.K_rr[tK.rAidx[sa[1][1]],tK.rBidx[sb[1][1]]]-\
    tK.K_lr[tK.lAidx[sa[1][0]],tK.rBidx[sb[1][1]]]-\
    tK.K_rl[tK.rAidx[sa[1][1]],tK.lBidx[sb[1][0]]])**2
    k_o=(tK.K_ll[tK.lAidx[sa[1][0]],tK.lBidx[sb[1][0]]]+\
    tK.K_rr[tK.rAidx[sa[1][1]],tK.rBidx[sb[1][1]]]+\
    tK.K_lr[tK.lAidx[sa[1][0]],tK.rBidx[sb[1][1]]]+\
    tK.K_rl[tK.rAidx[sa[1][1]],tK.lBidx[sb[1][0]]])       
    k_tppk=(tK.K_ll[tK.lAidx[sa[1][0]],tK.lBidx[sb[1][0]]]*\
    tK.K_rr[tK.rAidx[sa[1][1]],tK.rBidx[sb[1][1]]]+\
    tK.K_lr[tK.lAidx[sa[1][0]],tK.rBidx[sb[1][1]]]*\
    tK.K_rl[tK.rAidx[sa[1][1]],tK.lBidx[sb[1][0]]])
    kv=k_tppk+k_mlpk+k_o
    return kv

def getKernelData(K,S,dK=None,E=None):
    """
    Return the myKernelData object (convert to pyml style kernel)
    """
    #pdb.set_trace()
    mK=myKernelData(np.double(K),patternID=zip(*zip(*S)[:2]),labels=list(zip(*S)[-1]))        
    return mK
        
class complexKernel:
    """
    Evaluate the kernel values between 2 complexes A & B as follows
    (let l and r represent the ligand and receptor within each complex)
        K_ll=K(lA,lB)
        K_lr=K(lA,rB)
        K_rl=K(rA,lB)
        K_rr=K(rA,rB)
    It evaluates these kernels only over examples passed in E
    lAidx,rAidx,lBidx,rBidx are dictionary objects which when given the 
    residue index will return the index of that residue in the respective kernel
    Note: the inner kernel is computed using protKernel
    """
    def __init__(self,cidA,cidB,pdbpklpath,E,aAll=0,bAll=0):
        self.cidA=cidA
        self.cidB=cidB
        lAname=os.path.join(pdbpklpath,cidA+'_l_u.pdb.pkl')
        rAname=os.path.join(pdbpklpath,cidA+'_r_u.pdb.pkl')
        lBname=os.path.join(pdbpklpath,cidB+'_l_u.pdb.pkl')
        rBname=os.path.join(pdbpklpath,cidB+'_r_u.pdb.pkl')        
        lA=myPDB.loader(lAname)
        rA=myPDB.loader(rAname)
        lB=myPDB.loader(lBname)
        rB=myPDB.loader(rBname)
        #pdb.set_trace()
        if aAll:
            lAidx=range(len(lA.R))
            rAidx=range(len(rA.R))
        else:
            lAidx_p,rAidx_p=tuple(map(np.unique,zip(*E.Pex[cidA][0])))
            lAidx_n,rAidx_n=tuple(map(np.unique,zip(*E.Nex[cidA])))
            lAidx=list(np.unique(np.concatenate((lAidx_p,lAidx_n))))
            rAidx=list(np.unique(np.concatenate((rAidx_p,rAidx_n))))
        if bAll:            
            lBidx=range(len(lB.R))
            rBidx=range(len(rB.R))
        else:
            lBidx_p,rBidx_p=tuple(map(np.unique,zip(*E.Pex[cidB][0])))
            lBidx_n,rBidx_n=tuple(map(np.unique,zip(*E.Nex[cidB])))
            lBidx=list(np.unique(np.concatenate((lBidx_p,lBidx_n))))
            rBidx=list(np.unique(np.concatenate((rBidx_p,rBidx_n))))
        self.K_ll=protKernel(lA,lB,lAidx,lBidx).K
        self.K_lr=protKernel(lA,rB,lAidx,rBidx).K
        self.K_rl=protKernel(rA,lB,rAidx,lBidx).K
        self.K_rr=protKernel(rA,rB,rAidx,rBidx).K
        print "#NANS in kernels: ",np.sum(np.isnan(self.K_ll)),np.sum(np.isnan(self.K_lr)),np.sum(np.isnan(self.K_rl)),np.sum(np.isnan(self.K_rr))
        self.lAidx=dict(zip(lAidx,range(len(lAidx))))#lAidx;#
        self.rAidx=dict(zip(rAidx,range(len(rAidx))))#rAidx;#
        self.lBidx=dict(zip(lBidx,range(len(lBidx))))#lBidx;#
        self.rBidx=dict(zip(rBidx,range(len(rBidx))))#rBidx;#
    def __str__(self):
        return "complexKernel Object over complexes :"+self.cidA+','+self.cidB
        
def getValidE(E0,pdbpklpath):
    """
    
    """
    E=deepcopy(E0)
    K=set(E.Pex.keys())
    import os
    import glob
    s=set()
    for k in K:
        f=os.path.join(pdbpklpath,k+"*.pdb.pkl")
        if len(glob.glob(f))==4:
            s.add(k)
    tlist=list(s)
    E.Pex=copy_dict(E.Pex,*tlist)
    E.Nex=copy_dict(E.Nex,*tlist)
    return E
    
if __name__=="__main__":    
    pdbpklpath='./DBD4PDBPKL/' # '/s/chopin/b/grad/minhas/PDBPKLP3'#'/s/chopin/b/grad/minhas/Desktop/PDBPKLA' #
    exfname=pdbpklpath+'/E_6.0.lbl.pkl' #    '/media/sf_Desktop/PIANO/PDBPKL/E_6.0.lbl.pkl'# 
    ofname='dbdKernel.dbk.pkl'
    
    E=getExamplesDBD.loader(exfname)
    #pdb.set_trace()
    #skeys=E.Pex.keys()[0:2]
    #E.Pex=copy_dict(E.Pex,*skeys)
    #E.Nex=copy_dict(E.Nex,*skeys)
    tlist=['1A2K','1ACB','2A5T', '3CPH']
    E.Pex=copy_dict(E.Pex,*tlist)
    E.Nex=copy_dict(E.Nex,*tlist)
    

    dK=dbdKernel(pdbpklpath,E,ofname)
    
