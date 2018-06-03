# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:57:06 2013
Code to analyze the distance deviation for complexes
Distance deviation is the difference in the average distance of the top N
examples and randomly selected selected N examples that do not interact.
@author: root
"""
import random
from BISEPutils import *
from myPDB import *
from analyzePredFile import readFile,sortScores
from scipy.stats import nanmean, nanmedian
import traceback,os
import myPickle as Pickle
from PyML.evaluators import roc
import scipy
#def getDistMean(lD,r,top=True,N=10):
#    N=np.min((len(r),N))
#    if top:
#        idx=np.unique(r[:N])
#        N=len(idx)
#    else:
#        S=np.triu(lD)
#        S[S==0]=np.nan
#        return nanmean(nanmean(S))
##        idx=np.array(random.sample(r, np.min((len(r),3*N))))
#    d=[]
#    
#    for i in range(N):
#        for j in range(i+1,N):
#            d.append(lD[idx[i],idx[j]])
#    return np.mean(d)
def getDistMean(lD,r,top=True,M=40,N=40):
    """
    Compute the average distance between top M examples and remaining N examples
    
    """
    N=np.min((len(r),N))
    M=np.min((len(r),M))
    if top:
        idx=np.unique(r[:M])   
    else:
        tidx=np.unique(r[:M])
        idx=list(set(np.array(random.sample(r, np.min((len(r),N))))).difference(tidx))
    N=len(idx)
    print 'N=',N
    d=[]    
    for i in range(N):
        for j in range(i+1,N):
            d.append(lD[idx[i],idx[j]])
    #pdb.set_trace()            
    return np.mean(d)    
def computeDistMeansForComplex(cid,N,pdbpklpath,pppath):
    """
    code for getting distance and auc information
    """
    L=myPDB.loader(os.path.join(pdbpklpath,cid+'_l_u.pdb.pkl'))
    R=myPDB.loader(os.path.join(pdbpklpath,cid+'_r_u.pdb.pkl'))
    if type(pppath)==type(''):
        (pauc,Mv,Ml,lseq,rseq,lauc,rauc)=readFile(pppath+cid+'.pairpred.txt',usePDBidx=False)
    else:
        (pauc,Mv,Ml,lseq,rseq,lauc,rauc)=pppath
    lauc=None
    rauc=None
    try:
        (_,_,lauc)=roc.roc(list(np.array(lauc.values())[:,0]),list(np.array(lauc.values())[:,1]))
        (_,_,rauc)=roc.roc(list(np.array(rauc.values())[:,0]),list(np.array(rauc.values())[:,1]))
    except:
        pass
    Mlx=np.random.random(Ml.shape)
    Mlx[Ml<0]=-1
    (r,c,v)=sortScores(Mlx)
    
    lD=getDistMat(getCoords(L.R))
    rD=getDistMat(getCoords(R.R))    
    #pdb.set_trace()
    M=20
    return pauc,lauc,rauc,getDistMean(lD,r,top=True,M=M), getDistMean(lD,r,top=False,M=M), getDistMean(rD,c,top=True,M=M), getDistMean(rD,c,top=False,M=M)

def getRMSD():
    """
    get a dictionary object giving the rmsd, interface asa and category from the ./rmsd.txt file
    """
    frmsd='./rmsd.txt'
    drmsd={}
    for ln in open(frmsd):
        if ln[0]=='#':
            ln.split()
            categ=ln[1:].strip()
            continue
        lns=ln.split()
        cid=lns[0][:4]
        rmsd=float(lns[1])
        dasa=float(lns[2])
        if cid in drmsd:
            print cid,'Already there'
        drmsd[cid]=(rmsd,dasa,categ)
    return drmsd
def parallelRun(N,pdbpklpath,pppath,ofname=None,comm=None,myid=0,nprocs=1):
    """
    Wrapper for running computeDistMeansForComplex in parallel
    """
    cids=incids     
    csize=int(np.ceil(len(cids)/float(nprocs)))
    gclist=list(chunks(cids,csize))
    
    mycids=gclist[myid]
    A={}
    for cid in mycids:
        print "Processing",cid        
        try:
            A[cid]=computeDistMeansForComplex(cid,N,pdbpklpath,pppath)
            
        except Exception as ee:
            print "Error processing", cid,ee, traceback.format_exc()
            continue        
    Ascores=None
    if(myid!=0):
        comm.send(A, dest=0)        
    if(myid==0):
        gcK=[A]
        for p in range(1,nprocs):
            gcK.append(comm.recv(source=p))            
        Ascores=mergeDicts(gcK)
        if ofname is not None:            
            Pickle.dump(ofname,Ascores)
            print "Saved scores file:",ofname
    return Ascores       

if __name__=="__main__":
    ofname='../ppdists_esr_rand_true.dst.pkl'#_lbl
    pltt=True
    f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN',  '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
    incids=f3+f4
    #incids=incids[:1]
    
    if not(os.path.isfile(ofname) and os.access(ofname, os.R_OK)):
        bdir='../DBD4N/' # '/s/chopin/b/grad/minhas/Desktop/DBD4N/' # 
        pdbpklpath=bdir+'/PDBPKL4'
        pppath='../DBD4_ESR_prop/'

        N=20
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
        
        A=parallelRun(N,pdbpklpath,pppath,ofname=ofname,comm=comm,myid=myid,nprocs=nprocs)
    else:
        A=Pickle.load(ofname)
    if A is not None:        
        rmsd=getRMSD()    
        sels=['ALL','RB','MED','HARD'] #
        #sels=['ALL']
        mrks={'ALL':'o','RB':'o','MED':'s','HARD':'^'};
        clrs={'ALL':'b','RB':'g','MED':'b','HARD':'k'};    
        lbns={'ALL':'','RB':'Rigid Body','MED':'Medium','HARD':'Hard'};    
        for sel in sels:        
            
            if sel!='ALL':
                incids_sel=[cid for cid in rmsd.keys() if rmsd[cid][2]==sel]
            else:
                incids_sel=incids
            lbn=lbns[sel]
            mrk=mrks[sel]
            clr=clrs[sel]        
            As={}
            bsa={}
            for cid in incids_sel:
                if cid in A and cid in rmsd:
                    As[cid]=A[cid]
                    bsa[cid]=rmsd[cid][0]
            bsax=np.array(bsa.values())                    
            vA=np.array(As.values())
            auc=vA[:,0]
            lauc=vA[:,1]
            rauc=vA[:,2]
            pl=vA[:,3]
            dl=vA[:,4]
            pr=vA[:,5]
            dr=vA[:,6]
            Dl=pl-dl
            Dr=pr-dr
            print "Type of complexes:",sel
            print "Number of complexes",len(auc)
            
            
            xauc=np.array(list(lauc)+list(rauc))
            print "Mean AUC",nanmean(auc), nanmean(list(lauc)+list(rauc))
            xd=np.array(list(Dl)+list(Dr))
            print "Mean deviation",nanmean(xd)
            cc,pv=scipy.stats.pearsonr(xauc[~np.isnan(xd)],xd[~np.isnan(xd)])
            
            print "Correlation for Protein (dev,auc)",cc,pv
            mv=np.min((Dl,Dr),axis=0)
            cc,pv=scipy.stats.pearsonr(mv[~np.isnan(mv)],auc[~np.isnan(mv)])
            print "Correlation for complex (mindev,auc)",cc,pv
            mv=np.max((Dl,Dr),axis=0)
            cc,pv=scipy.stats.pearsonr(mv[~np.isnan(mv)],auc[~np.isnan(mv)])
            print "Correlation for complex (maxdev,auc)",cc,pv
            mv=np.mean((Dl,Dr),axis=0)
            cc,pv=scipy.stats.pearsonr(mv[~np.isnan(mv)],auc[~np.isnan(mv)])
            print "Correlation for complex (avgdev,auc)",cc,pv
        #    d=np.mean((vA[:,3]-vA[:,4],vA[:,5]-vA[:,6]),axis=0)
            l1=np.array(list(pr[~np.isnan(pr)])+list(pl[~np.isnan(pl)]))
            l2=np.array(list(dr[~np.isnan(pr)])+list(dl[~np.isnan(pl)]))
            idxx=~(np.isnan(l1)+np.isnan(l2))
            l1=l1[idxx]
            l2=l2[idxx]
            print "Wilcoxon Signed Rank Sum Test:",scipy.stats.wilcoxon(l1,l2)#scipy.stats.ttest_rel(l1,l2)
            
            print "WW: ",sel,nanmean(list(pl)+list(pr)),nanmean(list(dl)+list(dr))
           
#import myPickle
#rmsd0=myPickle.load(rmsdfname)
#Rx=[]
#for cid in As.keys():
#    try:                        
#        Rx.append([rmsd0[cid+'_l'][0],As[cid][3]-As[cid][4]])
#        Rx.append([rmsd0[cid+'_r'][0],As[cid][5]-As[cid][6]])
#    except:
#        continue
#Rx=np.array(Rx)    
#print spearmanr(Rx[:,0],Rx[:,1]),pearsonr(Rx[:,0],Rx[:,1])
                    
            if pltt and sel!='ALL':  
                plt.figure(0);plt.scatter(list(Dl)+list(Dr),list(lauc)+list(rauc),marker=mrk,c=clr,label=lbn); 
                plt.figure(1); plt.scatter(mv,auc,marker=mrk,c=clr,label=lbn);
                plt.figure(2);plt.scatter(list(pl)+list(pr),list(dl)+list(dr),marker=mrk,c=clr,label=lbn);
        if pltt:
            plt.figure(1);plt.legend(loc=0);plt.xlabel('Avg. Distance Deviation in complex');plt.ylabel('AUC for complex');plt.grid();#plt.savefig('../figures/fig_A11.png',  format='png', dpi=1200)

            plt.figure(0);plt.legend(loc=0);plt.xlabel('Distance Deviation');plt.ylabel('AUC for protein');plt.grid();#plt.savefig('../figures/fig_A10.png',  format='png', dpi=1200)

            plt.figure(2);plt.plot([0, 50],[0, 50],'r',linewidth=2);plt.legend(loc=0);plt.xlabel('Avg. distance between interacting residues');plt.ylabel('Avg. distance between randomly selected residues');plt.title('Protein level spatial clustering analysis');plt.grid();#plt.savefig('../figures/fig_A9.png',  format='png', dpi=1200)

            plt.show()