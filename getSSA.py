# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:09:37 2013
MI Calculation. For a single MSA file. For a pair of MSA files. For DBD complexes.
usage:  python getSSA.py [runs in parallel mode]
        python getSSA.py cid [runs for a single complex]
@author: root
"""
from myPDB import *
from getExamplesDBD import *
from hhblits import runHHblits
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
# import pdb
import sys
import numpy as np
import myPickle as cPickle

    
def getOX(s):    
    """
    get the organism code from the fasta header/description
    """
    try:
        s=s[s.find('OX=')+3:]
        s=s[:s.find(' ')]
    except:
        s=None
    return s
def getSDict(fr):
    """
    Given the MSA file, it creates a dictionary based on the species
    Input:
        fr: name of input fasta msa file
    Output:
        D: Dictionary with keys being the species code. The first entry is 
            given the code 'first'. If multiple alignments from a single 
            species are present in the input, the ones with the fewest gaps
            is retained.
    """
    handle = open(fr, "rU")
    first=True
    D={}
    for record in SeqIO.parse(handle, "fasta"):
        if first:
            spid='first'
            first=False
        else:
            spid=getOX(record.description)
            if spid is None:
                continue
        #pdb.set_trace()
        if (spid not in D) or (D[spid].seq.count('-')+D[spid].seq.count('.'))>(record.seq.count('-')+record.seq.count('.')): #Since the original records are sorted in increasing E-values, doing this retains the best match within a species
            D[spid]=record        
    handle.close()
    return D

def getCombinedMSA(Dl,Dr,fname=None):
    """
    Returns the combined MSA. The concatenated MSA is obtained in a species
    specific fashion (based on OX codes in the MSA)
    Input:
        Dl: Dictionary object obtained using getDict from the ligand MSA 
        Dr: Dictionary object obtained using getDict from the receptor MSA
        fname: (optional) When not None, it must be a string to which the combined MSA is saved in fasta format
    Return:
        Rlist: list of SeqRecord objects contained the combined alignment
        ll: length of the ligand sequence
    """
    ckeys=['first']
    ckeys.extend(list((set(Dr.keys()).intersection(Dl.keys()).difference(['first']))))
    rlist = []
    ll=0
    for k in ckeys:        
        desc=Dl[k].description+'##'+Dr[k].description
        if k=='first':
            desc='##'
            tt=Dl[k].seq.tostring()
            ll=len(tt)-(tt.count('-')+tt.count('.'))
        r=SeqRecord(Dl[k].seq+Dr[k].seq,id=Dl[k].id,description=desc)        
        rlist.append(r)
    if fname is not None:
        output_handle = open(fname, "w")
        SeqIO.write(rlist, output_handle, "fasta")
        output_handle.close()
    return rlist,ll

def calcMI(rlist,ll=0,fseq=0,x=0.7,AA='ACDEFGHIKLMNPQRSTVWY-'):
    """
    Computes mutual information given a MSA
    rlist: list of sequence records read from fastA file
    ll: length of first sequence in a combined MSA (when 0, then only 1 sequence is considered)
    fseq: focus sequence (sequence of interest, default -- the first in rlist)
    x: Threshold
    AA: amino acids
    Returns:
        Mi: Mututal information
        Fij: Represents Fij(A,B) -- the probability of having A@i and B@j
        Fi: Represents Fi(A) -- the probability of having A@i
        Cij: Represents Fij(A,B)=Fij(A,B)-Fi(A)*Fj()    
    """    
    ml=[]    
    for m in rlist:
        ml.append(list(m.seq.tostring()))
    ma=np.array(ml)
    ma=ma[:,(ma[fseq,:]!='-')*(ma[fseq,:]!='.')]   
    dAA=dict(zip(AA,range(len(AA))))
    q=len(AA)    
    M,L=ma.shape
    D=np.zeros((M,M))
    
    for i in range(M):
        D[i,i]=L
        for j in range(i+1,M):          
            D[i,j]=np.sum(ma[i,:]==ma[j,:])
            D[j,i]=D[i,j]            
    W=1.0/np.sum((D-x*L)>=0,axis=0)
    Meff=np.sum(W)
    
    Fi=np.zeros((L,q))
    for m in range(M):
        for i in range(L):
            try:
                Fi[i,dAA[ma[m,i]]]+=W[m]
            except KeyError: #if the character is unknown
                continue   
    
    Fij=np.zeros((ll,L-ll,q,q))    
    for m in range(M):
        for i in range(ll):            
            for j,jll in enumerate(range(ll,L)):
                try:
                    Fij[i,j,dAA[ma[m,i]],dAA[ma[m,jll]]]+=W[m]
                except KeyError: #if the character is unknown
                    continue
    Fij=Fij/Meff
    Fi=Fi/Meff
    Mi=np.zeros((Fij.shape[0],Fij.shape[1]))
    C=np.zeros((Fij.shape[0],Fij.shape[1],q,q))
    for i in range(Fij.shape[0]):
        for j in range(Fij.shape[1]):
            num=Fij[i,j]
            den=np.dot(np.atleast_2d(Fi[i]).T,np.atleast_2d(Fi[j+ll]))
            C[i,j]=num-den
            pidx=(num>0)*(den>0)
            if np.any(pidx):              
                Mi[i,j]=np.sum(num[pidx]*np.log(num[pidx]/den[pidx]))
    return (Mi,Fij,C,Fi,Meff)
    
def saveMIFile(M,ofile):
    """
    save the MI variables in a zipped pickle binary file
    """
    cPickle.dump(ofile, M)
def loadMIFile(ofile):
    return cPickle.load(ofile)
    
def calcMIForComplex(pdbpklpath,cid,odir='./',saveFile=False):
    """
    given a complex 'cid', it computes the MSA and MI
    Input:
        pdbpklpath: path where the feature files are kept
        cid: either strng complex id or tuple of (L,R) myPDB objects for that cid
        odir: (optional) where the output file will be saved. If this location already as
            <cid>.ec file, the no calculation is done and the MI is read from the file. If 
            no such file is found, then this function looks for the MSA files for the
            component proteins of cid and if those files are found then they are used, otherwise, 
            they are computed using hhblits.
        sabeFile: (optional) whether to saved the .ec file or not
    Output: 
        tuple: (mi,Fij,Cij,Fi,meff) where:
            mi is a numpy array of dimensions (lxr) where l and r are the lengths of L.R and R.R in the cid
                and it contains the mi values between each pairs. Locations that do not have valid amino acids
                have nans.
            Fij is a numpy array of dimension (lxrx21x21) and this represents Fij(A,B), i.e., the probability
                of observing A@i and B@j.
            Fi is a numpy array of dimensions ((len(L.seq)+len(R.seq))x21) representing Fi(A), i.e., the probability
                of observing A@i
            Cij(A,B) = Fij(A,B)-Fi(A)*Fj(B)
            meff: Number of effective sequences
        This tuple is saved in the file when saveFile is set to true
    """
    # Write the ligand
    if type(cid)==type(''):
        L=myPDB.loader(os.path.join(pdbpklpath,cid+'_l_u.pdb.pkl'))
        R=myPDB.loader(os.path.join(pdbpklpath,cid+'_r_u.pdb.pkl'))
    elif type(cid)==type(()):
        L=cid[0]
        R=cid[1]
        cid=L.name[:4]
    else:
        raise Exception("Unhandled complex input.")
    #pdb.set_trace()
    lfasta=os.path.join(odir,cid+'_l_u.seq')
    rfasta=os.path.join(odir,cid+'_r_u.seq')
    lafasta=lfasta+'.fasta'
    rafasta=rfasta+'.fasta'
    ecfile=os.path.join(odir,cid+'.mi.pkl')
    if os.path.exists(ecfile) and os.path.isfile(ecfile):#If the EC file already exists
        print "Using existing EC file", cid, ecfile
        return loadMIFile(ecfile) 
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
        M=calcMI(*getCombinedMSA(getSDict(lafasta),getSDict(rafasta)))
        #The code below sets the MI, Fij and C values of all the locations not in the sequence to nans
        nL=len(L.R)
        nR=len(R.R)
        q=M[1].shape[-1]        
        Mc=np.zeros((nL,nR)),np.zeros((nL,nR,q,q)),np.zeros((nL,nR,q,q)),M[-2],M[-1]
       
        for i in range(len(Mc)-2):
            Mc[i].fill(np.nan)
            for j,jr in enumerate(L.S2Ri):
                for k,kr in enumerate(R.S2Ri):
                    Mc[i][jr,kr]=M[i][j,k]
        
        
        if saveFile:
            saveMIFile(Mc,ecfile)
        return Mc
def parallelRun(pdbpklpath,exfname,odir,incids,comm=None,myid=0,nprocs=1):
    """
    For running the MSA evaluation in parallel to save the EC files
    Input:
        pdbpklpath: Where the feature files are kept
        exfname: path to the example file
        odir: ouput directory
        incids: list of cids of complexes on which evaluation is to be performed
        (comm, myid,nprocs): mpi4py parameters
    """
    from dbdKernels import getValidE
    E=getValidE(getExamplesDBD.loader(exfname),pdbpklpath,incids)
    cids=E.Pex.keys()
    csize=int(np.ceil(len(cids)/float(nprocs)))
    gclist=list(chunks(cids,csize))
    mycids=gclist[myid]
    for cid in mycids:
        try:
            calcMIForComplex(pdbpklpath,cid,odir,True)
        except:
            print "Error processing", cid
            continue
if __name__=="__main__":    
    if len(sys.argv)>1:
        cid=sys.argv[1]
        pdbpklpath='./DBD4N/PDBPKL4'#'/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'#        
        odir='./ecstats_rand2'
        calcMIForComplex(pdbpklpath,cid,odir,True)        
    else:
        pdbpklpath='/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'#'./DBD4N/PDBPKL4'#
        exfname=pdbpklpath+'/E_125PN_15_35_50.lbl.pkl'
        odir='./DBD4MSA'
        f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
        f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
        incids=f3+f4  
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            myid = comm.Get_rank()
            nprocs = comm.Get_size()
        except ImportError:
            print "Failure importing MPI4py: Not using MPI parallelization."
            comm=None
            myid=0
            nprocs=1        
        
        parallelRun(pdbpklpath,exfname,odir,incids,comm=comm,myid=myid,nprocs=nprocs)
    

    