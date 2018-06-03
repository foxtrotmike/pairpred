# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 12:03:18 2013
Leave one complex out cross validation using the transductive pairwise constraint svm
@author: root
"""
from testComplex_PWCC import *
if __name__=="__main__":
    cargs={'nCycles': 3000000,'Lambda':3e-6,'w':0.002,'verbose':2,'nConstraints':30000,'iCycles':500000, 'pUpdate':0.0,'eUpdate':3} #default arguments for the classification
    pdbpklpath='../../DBD4CSPKL/PKL/'
    pwKfname='dbdKernel.dbk.pkl'
    odir='./'
    evalROC=True #Whether we want to calculate the ROC or not 
    selectAll=True  #Whether to use all examples in the complex or not
    bsize=3000
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    nprocs = comm.Get_size()
    cids=['2OOB',  '1PPE', '1J2J', '1GL1', '1SYX', '1Z0K', '1AY7', '1FFW', '3SGQ', '1S1Q', '1FLE', '7CEI', '2IDO',  '4CPA', '2UUY', '1R6Q', '1D6R', '1OC0', '1CGI', '1R0R', '1EAW',  '1XD3', '1LFD', '2I25', '1CLV', '1H9D', '1ACB', '2SNI', '3D5S', '1Z5Y', '2HRK', '2ABZ', '1UDI', '1PXV', '2J0T']#E.Pex.keys()[20:40]    '1EFN','1KTZ','1GCQ',
    N=float(len(cids))
    bindx=int(np.floor(myid*N/nprocs))
    eindx=int(np.floor((myid+1)*N/nprocs))
    mycids=cids[bindx:eindx]    
    print 'I am process:',myid,' and complexes being processed by me are:',len(mycids),mycids    
    for (i,cid) in enumerate(mycids):        
        fname=odir+('#'.join(cid))+'.pairpred.txt'
        print 'Currently Processing: ',i,cid
        try:
            testComplex(cid,pwKfname,pdbpklpath,evalROC=evalROC,selectAll=selectAll,ofile=fname,bsize=bsize,**cargs)
        except Exception as e:
            print '###PROCESSSING FAILED FOR ',i,cid,e
        
                