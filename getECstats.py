# -*- coding: utf-8 -*-
"""
Created on Fri May  3 00:00:18 2013
Run hhblits on all complexes in the list "incids" and compute MI and DI using
calcEC.m. It also calculates the AUC scores.
@author: root
"""

#load E
#for each complex in E, compute its DI and MI
#Save these as [cid,lrid,rrid,lbl,di,mi]
from calcMIDI import getECForComplex, getAUCs, toArray
from dbdKernels import getValidE
from getExamplesDBD import *
import cPickle
import glob
from getSSA import getSDict
#from DICOR import getECForProtein,getAUCs
from BISEPutils import getFileParts
if __name__=="__main__":    
    pdbpklpath='./DBD4N/PDBPKL4'#'/s/chopin/b/grad/minhas/Desktop/DBD4N/PDBPKL4'
    exfname=pdbpklpath+'/EP_6N.lbl.pkl'
    odir='./'
    ofname='ecstats_rand.pkl'
    f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
    incids=f3+f4 #None
    
    #incids=incids[:]
    #incids=['1SBB','1H9D','1AY7']
    #incids=[getFileParts(f)[1] for f in glob.glob(os.path.join(odir,'*.ec'))]
    incids=['2I25', '2MTA', '2NZ8', '2O3B', '2O8V', '2OOR', '2OZA', '2SNI', '3D5S', '3SGQ', 'BOYV', '1QFW', '1R0R', '1SBB', '1T6B', '1US7', '1WDW', '1WEJ', '1Z0K', '2BTF', '2CFH', '2H7V', '1IQD', '1JMO', '1JWH', '1K4C', '1KKL', '1KXP', '1M10', '1MLC', '1MQ8', '1N8O', '1NSN', '1OYV', '1FSK', '1H9D', '1HCF', '1HE1', '1HE8', '1HIA', '1I2M', '1I4D', '1IJK', '1ACB', '1AK4', '1AY7', '1BJ1', '1BKD', '1DQJ', '1E96', '1EER', '1F34', '1F6M', '1FC2', '1FCC']
    #incids=['2I25', '2MTA', '2NZ8', '2O3B', '2O8V', '2OOR', '2OZA', '2SNI', '3D5S', '3SGQ', 'BOYV', '1QFW', '1R0R', '1SBB']
    #incids=['1SBB']
    #incids=incids[:1]
    incids=['1AY7']
    evalROC=True
    try:
        from PyML.evaluators import roc        
    except ImportError:
        print "PyML not found, will not evaluate ROC."
        evalROC=False 
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        nprocs = comm.Get_size()
    except ImportError:
        print "Failure importing MPI4py: Not using MPI parallelization"
        comm=None
        myid=0
        nprocs=1
    
    E=getValidE(getExamplesDBD.loader(exfname),pdbpklpath,incids)
    cids=E.Pex.keys()
    csize=int(np.ceil(len(cids)/float(nprocs)))
    gclist=list(chunks(cids,csize))
    mycids=gclist[myid]
    myEdata=[]
    Lauc=[]
    Xauc=[]
    Xuauc=[]
    NCS=[]
    ncs=0
    for cid in mycids:
        (mx,dx,dlr),lstats,rstats=getECForComplex(pdbpklpath,cid,odir=odir)
        #ncs=len(set(getSDict(os.path.join(odir,cid+'_l_u.seq.fasta')).keys()).intersection(getSDict(os.path.join(odir,cid+'_r_u.seq.fasta')).keys()))
        NCS.append(ncs)
        mxa=toArray(mx,lstats[0].shape[0],rstats[0].shape[0])
        pdb.set_trace()
        dxa=toArray(dx,lstats[0].shape[0],rstats[0].shape[0])
        Xuauc.append(getAUCs(mxa,dxa,dlr,dthr=6.0))
        Lauc.append(getAUCs(*rstats))
        Lauc.append(getAUCs(*lstats))        
        dd=[]
        ld=[]
        md=[]
        for p in E.Pex[cid][0]:
            if p in mx:
                #myEdata.append([cid,p[0],p[1],1,mx[p],dx[p]])
                dd.append(dx[p])
                ld.append(+1)
                md.append(mx[p])
        for n in E.getNegEx(cid):
            if n in mx:
                #myEdata.append([cid,n[0],n[1],-1,mx[n],dx[n]])
                dd.append(dx[n])
                md.append(mx[n])
                ld.append(-1)
        (_,_,aa_di)=roc.roc(dd,ld)
        (_,_,aa_mi)=roc.roc(md,ld)
        Xauc.append([aa_mi,aa_di])
        print cid,ncs,lstats[0].shape[0],rstats[0].shape[0],Xauc[-1],Xuauc[-1],Lauc[-2:]
        #pdb.set_trace()
    if(myid!=0):
        comm.send(myEdata, dest=0)
    else:
        """
        for p in range(1,nprocs):
            myEdata.extend(comm.recv(source=p))
        output = open(ofname, 'wb')
        cPickle.dump(myEdata, output,-1)        
        output.close()     
        if evalROC:
            MV=[]
            DV=[]
            for cid in cids:
                dx=[di for (c,el,er,lbl,mi,di) in myEdata if c==cid]
                mx=[mi for (c,el,er,lbl,mi,di) in myEdata if c==cid]
                lx=[lbl for (c,el,er,lbl,mi,di) in myEdata if c==cid]
                MV.append((list(mx),list(lx)))
                DV.append((list(dx),list(lx)))
            (fp_mi,tp_mi,auc_mi)=roc.roc_VA(MV)
            (fp_di,tp_di,auc_di)=roc.roc_VA(DV)
            print 'AUC_MI =',auc_mi
            print 'AUC_DI =',auc_di
        """
        print 'Xauc =',np.mean(Xauc,axis=0)
        print 'Local AUC =',np.mean(Lauc,axis=0)
        print 'UAUC =',np.mean(Xuauc,axis=0)