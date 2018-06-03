# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 20:57:34 2013
Develop a dictionary of amino acid frequencies in NRPDB using SRF files
@author: root
"""
from BISEPutils import *
import os,sys,glob,tempfile,shutil
from mpi4py import MPI
import numpy as np
PROBIS_NRPDB_PATH='/s/parsons/w/nobackup/protfun/pairpred/NRPDB'#'../../NRPDB'#

def countAAinSRF(sfile):
    """
    Computes the number of times each amino acid appears in the SRF file
    """
    S={}
    with open(sfile,'r') as s:
        for l in s:
            if len(l)>1 and l[:2]=='A>':
                ls=l.split()[8:10]
                S[ls[1]]=ls[0]
    C=np.zeros(len(aa3idx))
    for s in S.values():
        try:
            C[aa3idx[s]]+=1
        except KeyError:
            continue
    return C

def getAACinNRPDB(srfdir,comm=None,myid=0,nprocs=1):
    flist=None
    if myid==0:
        flist=glob.glob(os.path.join(srfdir,'*.srf'))
    if comm is not None:
        flist=comm.bcast(flist,root=0)
    
    N=float(len(flist))
    bindx=int(np.floor(myid*N/nprocs))
    eindx=int(np.floor((myid+1)*N/nprocs))
    myflist=flist[bindx:eindx]
    C=np.zeros(len(aa3idx))
    for f in myflist:
        try:
            C=C+countAAinSRF(f)   
        except:
            continue
    if myid==0:
        for i in range(1,nprocs):            
            C+=comm.recv(source=i)            
        print C
    else:
        comm.send(C,dest=0)
    
    #comm.Barrier() 
    
        
if __name__=="__main__":
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    nprocs = comm.Get_size()    
    getAACinNRPDB(PROBIS_NRPDB_PATH,comm=comm,myid=myid,nprocs=nprocs)