# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:07:24 2013
Code for comparing the performance of LOOCV over all examples with the one 
over selected examples
@author: root
"""
from BISEPutils import getFileParts
from analyzePredFile import *
from dbdscrpp3 import getAUC
import glob
from PyML.evaluators import roc   
from scipy.spatial.distance import cdist
# import pdb  
from scipy import stats
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
    pp=np.sum(Mlx==1) # total number of positives
    nn=len(Mlx)-pp #total number of negatives
    return (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx)
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
    LV=[]   
    TPS=[]
    DNTP=[]
    auconly=False # whether to calculate the avg. auc of the complexes or do more
    (auc,(fp,tp),(A,Rx,D,L,cids,r,dkey))=getAUC('./Results/result_dbd4n_tppk.res.pkl')
    dAo=dict(zip(cids,A)) #AUCs from the training data set (CV) only
    loopath='./DBD4_ESR_prop/'#C:\Users\Afsar\Desktop\pairpred\sequence only\DBD3LOOCVSEQ
    fs=glob.glob(loopath+'*.pairpred.txt')
    dA={}
    dsA={}
    #daoo={}
    for i,ifile in enumerate(fs):
        cid=getFileParts(getFileParts(ifile)[1])[1]
        print 'cid =',cid,100*float(i+1)/len(fs),'% done'
        #
        if auconly:
            auc=readFile(ifile,auconly=True)
            #aucoo=readFile('./DBD3LOOCV/'+getFileParts(ifile)[1]+getFileParts(ifile)[2],auconly=True)
            #daoo[cid]=aucoo
        else:
            (auc,ttp,fpi,dntp,la,ra,pp,nn,Mvx,Mlx)=computeNTP(ifile,top=200)            
            TPS.append([ttp,fpi,100.0*ttp/pp,pp,nn,pp+nn])
            DNTP.append(dntp)
            LV.append((list(Mvx),list(Mlx)))       
            dsA[cid]=(auc,la,ra)            
        dA[cid]=auc
        #print cid,auc,dAo[cid]
    print 'Number of complexes',len(dA)
    print 'Complex wise AUC = ',np.mean(dA.values()),'AUC for reduced set = ',np.mean([dAo[k] for k in dA.keys() if k in dAo.keys()])
    if not auconly:
        p12=map(list,zip(*dsA.values()));pa=p12[0];p1=p12[1];p2=p12[2];ps=p1;ps.extend(p2);
        print 'Complex Wise AUC =',np.mean(pa),'Protein Wise AUC =',np.mean(ps)  
        #ROC CURVE
        (fplv,tplv,auclv)=roc.roc_VA(LV)
        plt.figure();plt.plot(fplv,tplv);plt.xlabel('FP');plt.ylabel('TP');plt.grid();plt.title('ROC Curve: AUC =  %1.2f' % (auclv*100))
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
        plt.figure();plt.boxplot(tuple(XX),bootstrap=1000,positions=dthresh);plt.xlabel('Sequence Distance (D) from a TP'); plt.ylabel('Minimum rank of a prediction within distance D of a TP' );plt.title('Results of soft sequence distance threshold');plt.grid();plt.yticks(range(0,201,10));
        plt.show() 
       