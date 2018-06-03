# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 09:05:21 2013

@author: root
"""

from analyzeLOOCV_par_pwgen import *
import numpy as np
from postProcess import postProcessAvg
import glob,pdb,traceback
from BISEPutils import getFileParts 
import myPickle,re
pdbpklpath='../../DBD4CSPKL/PKL/'#'../../DBD4N/PDBPKL4'

#dirs=['../../DBD4_ESR_prop/','../../DBD4_NoESR_prop/','../../DBD4_SGDSVM/','../../DBD_SSVM_3E-6/','../../DBD_SPW_3E-6/']
#dirs=['../DBD4_ESR_prop/']#,'../SGD_DBD4/'
#dirs=['../DBD4_SGDSVM/','../DBD_SPW_3E-8/']
#cids=['2OOB', '1EFN', '1PPE', '1J2J', '1GL1', '1SYX', '1Z0K', '1AY7', '1FFW', '3SGQ', '1S1Q', '1FLE', '7CEI', '2IDO', '1KTZ', '4CPA', '2UUY', '1R6Q', '1D6R', '1OC0', '1CGI', '1R0R', '1EAW', '1GCQ', '1XD3', '1LFD', '2I25', '1CLV', '1H9D', '1ACB', '2SNI', '3D5S', '1Z5Y', '2HRK', '2ABZ', '1UDI', '1PXV', '2J0T']#E.Pex.keys()[20:40]    
#dirs=['../../SGD_DBD4/']
dirs=['../../DBD4S_SMO196/','../../DBD4S_SGD196/','../../DBD4S_SGD196_CENT/','../../DBD4S_SGDCENTPW71/']#['../../DBD4_NoESR_prop/']
# FIND THE COMMON SET 
cids=None
Re=(r"(\S+)\.",r"(\S+)\#(\S+)\.")
for d in dirs:
    cids_d=[]
    for f in glob.glob(d+'/*.pairpred.txt'):
        cx=re.match(Re[int('#' in f)],getFileParts(f)[1]).groups()
        if len(cx)==1:
            cx=cx[0]
        cids_d.append(cx)
    if cids is None:
        cids=set(cids_d)
    else:
        cids=cids.intersection(cids_d)
cids=list(cids)        
R={}
for cid in cids:
    rfpp=[]
    for d in dirs:    
        fname=d+('#'.join(cid))+'.pairpred.txt'
        try:            
            (auc0,Mv0,Mv,Ml,lseq,rseq,lrV0,lrV,rrV0,rrV)=postProcessAvg(cid,pdbpklpath,d) 
            nidx=~np.isnan(Mv) &  ~np.isnan(Ml)
            (_,_,auc)=roc.roc(list(Mv[nidx]),list(Ml[nidx]))
            rp=getTP_RFPP(Mv,Ml);
            rnp=getTP_RFPP(Mv0,Ml);
            pap=auc,getAUC4Protein(lrV)[0],getAUC4Protein(rrV)[0]
            panp=auc0,getAUC4Protein(lrV0)[0],getAUC4Protein(rrV0)[0]
            r=(pap+rp,panp+rnp)#+(100*np.mean(Ml[~np.isnan(Ml)]>0),Mv.shape[0],Mv.shape[1])
            #(auc,Mv,Ml,lseq,rseq,lrV,rrV)=readFile(fname)
        except Exception as e:
            print e
            print '-'*60
            print '###PROCESSSING FAILED FOR ',cid,e
            traceback.print_exc(file=sys.stdout)
            print '-'*60            
            r=np.nan
        rfpp.append(r)
    print cid,rfpp
    R[cid]=rfpp


# (AUC, AUCL, AUCR, NTP,RFP)_post,(AUC, AUCL, AUCR, NTP,RFP)_no_post,pp, percentage of positives, |L|,|R|
myPickle.dump('DBD4_SGD_CENTPW71.res.pkl',R)    
#V=np.array([r for r in R.values() if ~np.any(np.isnan(r))])    
#import scipy.stats
#scipy.stats.wilcoxon(V[:,0]-V[:,1])

#cc=['1KTZ',  '2OOB']

#Ro={};
#for k in R:
#    if k not in cc:
#        Ro[k]=[]
#        for m in range(len(R[k])):
#            Ro[k].extend(R[k][m][0]+R[k][m][1])
#        
#mV=np.mean(Ro.values(),axis=0)        
#print mV
#
#with open('eggs38.csv', 'wb') as csvfile:
#        spamwriter = csv.writer(csvfile, delimiter=',')
#        for cid in Ro:
#            spamwriter.writerow([cid]+list(Ro[cid]))    