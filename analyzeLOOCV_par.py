# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:07:24 2013
Code for comparing the performance of LOOCV over all examples with the one 
over selected examples
@author: root
"""
from BISEPutils import getFileParts,chunks,mergeDicts
from analyzePredFile import *
from dbdscrpp3 import getAUC
import glob
from PyML.evaluators import roc   
from scipy.spatial.distance import cdist
# import pdb  
from scipy import stats
from postProcess import postProcessAvg
import traceback
import myPickle as mkl
#import yard
def getAUC4Protein(lrV):
    vl=map(list, zip(*lrV.values()));vv=vl[0];ll=vl[1]    
    (_,_,a)=roc.roc(vv,ll)
    vv=np.array(vv)
    ll=np.array(ll)
    return (a,vv,ll)
def findNTPinTop(Mvx,Mlx,Mvshape,top):
    # returns : 
    #   ttp: Number of true positives in top
    #   fpi: index of first true positive
    #   dntp: distance to the nearest positive for all 'top' examples
    #find the number of true positives in the top 'top'
    sidx=np.argsort(Mvx)[::-1]
    L=[si for si,i in enumerate(sidx[:top]) if Mlx[i]==1]
    dntp=[] #distance from the nearest true positive
    #find the rank of the first positive
    
    (rv,cv)=np.unravel_index(sidx[:top], Mvshape)
    (rl,cl)=np.unravel_index(np.nonzero(Mlx==1), Mvshape)
    rcl=np.array((rl.flatten(),cl.flatten()))
    rcv=np.array((rv.flatten(),cv.flatten()))
    D=cdist(rcl.T, rcv.T, 'chebyshev')#cityblock
    di=np.argmin(D,axis=0)
    dntp=np.array([D[dix,i] for i,dix in enumerate(di)])    
    if len(L):
        fpi=L[0]        
    else:
        for si,i in enumerate(sidx[top:]):
            if Mlx[i]==1:
                break
        fpi=top+si
            
    ttp=len(L) #top true positives

    
    return (ttp,fpi,dntp)
    
def computeNTP(ifile,top=200):
    """
    Given a result file, it computes
    auc: The auc score
    ttp: Number of true positives in top
    fpi: index of the first true positive
    dntp: Distance to the nearest true positive for each top example
    la: auc of ligand
    ra: auc of receptor
    pp: number of positive examples
    nn: number of negative examples
    Mvx: flattened matrix of prediction scores
    Mlx: flattened matrix of labels
        
    """
    if type(ifile)==type(''):
        (auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(ifile,usePDBidx=False);#(auc,Mv,Ml,lseq,rseq,lrV,rrV)
    else: #expects tuple
        (auc,Mv,Ml,lseq,rseq,lrV,rrV)=ifile
    (la,lv,ll)=getAUC4Protein(lrV)
    (ra,rv,rl)=getAUC4Protein(rrV)
    Mvx=Mv.ravel()
    Mlx=Ml.ravel()
    nidx=~np.isnan(Mvx) &  ~np.isnan(Mlx)
    Mvx[~nidx]=-np.inf            
    (ttp,fpi,dntp)=findNTPinTop(Mvx,Mlx,Mv.shape,top=top)
    Mvx=Mvx[nidx]
    Mlx=Mlx[nidx]
    
    #yard.ROCCurve(yard.BinaryClassifierData(zip(Mvx,Mlx)))#PrecisionRecallCurve
    #zxA=yard.ROCCurve(yard.BinaryClassifierData(zip(Mvx,Mlx)))
    #zxR=yard.ROCCurve(yard.BinaryClassifierData(zip(rv,rl)))
    #zxL=yard.ROCCurve(yard.BinaryClassifierData(zip(lv,ll)))
    #pdb.set_trace()
    pp=np.sum(Mlx==1) # total number of positives
    nn=len(Mlx)-pp #total number of negatives
    #pdb.set_trace()
    return (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx,lv,ll,rv,rl)
def calcRFPP(FPI,DNTP,dthresh=[0,1,2,3,4],pctiles=[10,25,50,75,90]):
    """
    FPI: List of index of first true positive for multiple complexes
    DNTP: list of lists of distances of top examples for each complex from their nearest true positives
    
    """
    #DISTRIBUTION PLOT
     #percentiles chosen for the analysis
     # sequence distance threshold    
    XX=[]    
    XX.append(FPI)
    print  [stats.scoreatpercentile(FPI,p) for p in pctiles]
    for dx in dthresh[1:]:
        DD=[]
        for dn in DNTP:
            d=np.nonzero(dn<=dx)[0]
            if len(d):
                DD.append(d[0]+1)
            else:
                DD.append(200)
        XX.append(DD)
        print [stats.scoreatpercentile(DD,p) for p in pctiles]
    return XX
if __name__=='__main__':
    ofname='DBD3_ESR_prop_postprocess.mkl.roc'
    doplot=False
    postprocess=True
    auconly=False # whether to calculate the avg. auc of the complexes or do more
    loopath='../DBD4_ESR_prop/'#C:\Users\Afsar\Desktop\pairpred\sequence only\DBD3LOOCVSEQ
    bdir= '../DBD4N/' #'/s/chopin/b/grad/minhas/Desktop/DBD4N/' # 
    f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
    incids=f3 #['1SBB', '1JPS', '2HMI', '1GHQ', '1K74', '1D6R', '2SIC', '2HRK', '1EAW', '1VFB', '7CEI', '1I4D', '1H1V', '2PCC', '1FQ1', '1FQJ', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1GP2', '1FAK', '1NW9', '1GRN', '1A2K', '1JMO', '1RLB', '1DQJ', '2BTF', '2I25', '1BUH', '1BGX', '1EFN', '1DFJ', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1QFW', '2H7V', '1KAC', '1QA9', '1AHW', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1B6C', '1PPE', '2O8V', '1HIA', '2SNI', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IB1', '1GCQ', '1F51', '2B42', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IRA', '1E6J', '1Z5Y', '1EZU', '1E6E', '2C0L', '1TMQ', '1ZHI', '1M10', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '1Z0K', '1BJ1', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
    pdbpklpath=bdir+'/PDBPKL4'
    fs=glob.glob(loopath+'*.pairpred.txt')
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
    csize=int(np.ceil(len(fs)/float(nprocs)))
    gclist=list(chunks(fs,csize))    
    myfs=gclist[myid]
    LV=[]   
    TPS=[]
    DNTP=[]
    dsA={}
    LVP=[]
    for i,ifile in enumerate(myfs):        
        try:
            cid=getFileParts(getFileParts(ifile)[1])[1][:4]
            if incids is not None and cid not in incids:
                continue
            print "Processing",cid
            
            if auconly:
                auc=readFile(ifile,auconly=True)
                dsA[cid]=(auc,np.nan,np.nan)     
            else:                
                if postprocess:
                    pauc,Mvc0,Mvc,Mlc,lseq,rseq,lrV0,lrV,rrV0,rrV=postProcessAvg(cid,pdbpklpath,loopath)
                    ifile=(pauc,Mvc,Mlc,lseq,rseq,lrV,rrV)
                (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx,lv,ll,rv,rl)=computeNTP(ifile,top=200)    #lv,ll,rv,rl        
                TPS.append([ttp,fpi,100.0*ttp/pp,pp,nn,pp+nn])
                DNTP.append(dntp)
                LV.append((list(Mvx),list(Mlx))) 
                LVP.append((list(lv),list(ll))) 
                LVP.append((list(rv),list(rl))) 
                dsA[cid]=(auc,la,ra)            
        except Exception as ee:
            print "Error processing", cid,ee, traceback.format_exc()
            continue        
    if(myid!=0):
        comm.send((DNTP,dsA,LV,LVP,TPS), dest=0)        
    if(myid==0):
        dsAr=[dsA]
        for p in range(1,nprocs):
            (DNTP_p,dsA_p,LV_p,LVP_p,TPS_p)=comm.recv(source=p)
            dsAr.append(dsA_p)
            DNTP.extend(DNTP_p)
            LV.extend(LV_p)
            LVP.extend(LVP_p)
            TPS.extend(TPS_p)
        dsA=mergeDicts(dsAr)
        print 'Number of complexes',len(dsA)
        #print 'Complex wise AUC = ',np.mean(dA.values())
        p12=map(list,zip(*dsA.values()));pa=p12[0];p1=p12[1];p2=p12[2];ps=p1;ps.extend(p2);
        print 'Complex Wise AUC =',np.mean(pa),'Protein Wise AUC =',np.mean(ps)  
        if not auconly:
            (fplv,tplv,auclv)=roc.roc_VA(LV) 
            (fplvp,tplvp,auclvp)=roc.roc_VA(LVP) 
            mkl.save(ofname,((fplv,tplv,auclv),(fplvp,tplvp,auclvp)))
            print "AUC = ",auclv
            """        
                plt.hist(np.array(DNTP).flatten(),[0,1,2,3,4,5,6,1000],cumulative=True);plt.grid();plt.xlabel('sequence distance');plt.ylabel('counts');plt.title('Number of top 200 predictions vs. sequence distance from nearest true positive');plt.show()
                [np.sum(dn<2.0) for dn in DNTP]
                cids=[getFileParts(getFileParts(ifile)[1])[1] for ifile in fs]
                [dsA[cid] for cid in cids]
                [dAo[cid] for cid in cids]
            """ 
            #DISTRIBUTION PLOT
            dthresh=[0,1,2,3,4] # sequence distance threshold    
            XX=calcRFPP(np.array(TPS)[:,1]+1,DNTP,dthresh=dthresh)
            if doplot:
                plt.figure();plt.plot(fplv,tplv);plt.xlabel('FP');plt.ylabel('TP');plt.grid();plt.title('ROC Curve: AUC =  %1.2f' % (auclv*100))
                plt.figure();plt.boxplot(tuple(XX),bootstrap=1000,positions=dthresh);plt.xlabel('Sequence Distance (D) from a TP'); plt.ylabel('Minimum rank of a prediction within distance D of a TP' );plt.title('Results of soft sequence distance threshold');plt.grid();plt.yticks(range(0,201,10));
                plt.show() 
       