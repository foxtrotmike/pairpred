# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:07:24 2013
Generating the results AUC and RFPP from analyze1SVM
@author: root
"""

from BISEPutils import getFileParts
from getExamplesDBD import *
from analyzePredFile import *
from dbdscrpp3 import getAUC
import glob
from PyML.evaluators import roc   
from scipy.spatial.distance import cdist
# import pdb  
from scipy import stats
from scipy import sparse
def getAUC4Protein(lrV):
    vl=map(list, zip(*lrV.values()));vv=vl[0];ll=vl[1]    
    (_,_,a)=roc.roc(vv,ll)
    vv=np.array(vv)
    ll=np.array(ll)
    return (a,vv,ll)
def findNTPinTop(Mvx,Mlx,Mvshape,top):
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
def parseContLine(ln):
    # ['A', '#5', 'ASN:7', 'N', '::', 'B', '#5', 'HIS:6', 'H:', '0', '53.61']
    #   0   1       2       3     4     5   6       7       8    9      10
    lns=ln.split()
    lidx=lns[0]+lns[1]
    ridx=lns[5]+lns[6]
    lbl=int(lns[9])
    return (lidx,ridx,lbl)
def parse1SVM(cid,E,Asgl):
    (la,ra,lrV,rrV)=Asgl[cid]
    
    I=[]
    J=[]
    V=[]
    L=[]
    for xr in lrV.keys():
        for xc in rrV.keys():
            if (xr,xc) in E.Pex[cid][0]:
                l=+1.0
            else:
                l=-1.0
            I.append(xr)
            J.append(xc)
            v=lrV[xr][0]+rrV[xc][0]
            V.append(v)
            L.append(l)
    Mv=np.zeros((np.max(I)+1,np.max(J)+1))
    Mv.fill(np.nan)
    Ml=np.zeros(Mv.shape)    
    for idx in range(len(I)):
        Mv[I[idx],J[idx]]=V[idx]
        Ml[I[idx],J[idx]]=L[idx]
    return (Mv,Ml)
if __name__=='__main__':
    exfname='EP_6N.lbl.pkl'
    sglfile='result.sgl.pkl' 
    E=getExamplesDBD.loader(exfname)       
    Asgl=cPickle.load(open(sglfile, "rb" ) )

    LV=[]   
    TPS=[]
    DNTP=[]
    dA={}
    dsA={}
    #daoo={}

    
    for i,cid in enumerate(Asgl.keys()):
        
        print 'cid =',cid,100*float(i+1)/len(Asgl),'% done'
        (Mv,Ml)=parse1SVM(cid,E,Asgl)
        if ~np.any(Ml==1):
            continue
        #(auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(ifile,usePDBidx=False);#(auc,Mv,Ml,lseq,rseq,lrV,rrV)  
        
        #(la,lv,ll)=getAUC4Protein(lrV)
        #(ra,rv,rl)=getAUC4Protein(rrV)
        Mvx=Mv.ravel()
        Mlx=Ml.ravel()
        nidx=~np.isnan(Mvx) &  ~np.isnan(Mlx)
        Mvx[~nidx]=-np.inf            
        (ttp,fpi,dntp)=findNTPinTop(Mvx,Mlx,Mv.shape,top=500)
        Mvx=Mvx[nidx]
        Mlx=Mlx[nidx]
        pp=np.sum(Mlx==1) # total number of positives
        nn=len(Mlx)-pp #total number of negatives
        TPS.append([ttp,fpi,100.0*ttp/pp,pp,nn,pp+nn])
        DNTP.append(dntp)
        LV.append((list(Mvx),list(Mlx)))  
        
        #dsA[cid]=(auc,la,ra)
        #pdb.set_trace()    
        #dA[cid]=auc
        #print cid,auc,dAo[cid]
    #print 'Number of complexes',len(dA)
    #print 'Complex wise AUC = ',np.mean(dA.values()),'AUC for reduced set = ',np.mean([dAo[k] for k in dA.keys() if k in dAo.keys()])
    
    #p12=map(list,zip(*dsA.values()));pa=p12[0];p1=p12[1];p2=p12[2];ps=p1;ps.extend(p2);
    #print 'Complex Wise AUC =',np.mean(pa),'Protein Wise AUC =',np.mean(ps)  
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
    pctiles=[10,25,50,75,90] #percentiles chosen for the analysis
    dthresh=[0,1,2,3,4] # sequence distance threshold
    XX=[]    
    XX.append(list(np.array(TPS)[:,1]+1))
    print     [stats.scoreatpercentile(np.array(TPS)[:,1]+1,p) for p in pctiles]
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
    plt.figure();plt.boxplot(tuple(XX),bootstrap=1000,positions=dthresh);plt.xlabel('Sequence Distance (D) from a TP'); plt.ylabel('Minimum rank of a prediction within distance D of a TP' );plt.title('Results of soft sequence distance threshold');plt.grid();plt.yticks(range(0,201,10));
    plt.show() 
       