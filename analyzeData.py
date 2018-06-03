# -*- coding: utf-8 -*-
"""
Created on Tue Oct 09 21:40:36 2012
Check different properties of the data
@author: Afsar
"""

from myPDB import *
# import pdb
import glob
from getExamplesDBD import *

pklpath='../DBD4N/PDBPKL4/'
#fnames=glob.glob(pklpath+'*_u.pdb.pkl')
E=getExamplesDBD.loader(pklpath+'E_125PN_15_35_50.lbl.pkl')
cids=E.Pex.keys()
cids=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN',  '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']

A=[]
PP=[]
S=0
Sr=0
thr=0.05
for cid in cids:
    print cid
    try:
        L=myPDB.loader(pklpath+cid+'_l_u.pdb.pkl')
        R=myPDB.loader(pklpath+cid+'_r_u.pdb.pkl')
        lidx=set([])
        ridx=set([])
        P=0
        for p in E.Pex[cid][0]:
            lidx.add(p[0])
            ridx.add(p[1])
            if L.rASA[p[0]]>thr and R.rASA[p[1]]>thr:
                P=P+1
        PP.append((P,len(E.Pex[cid][0])))
        lASA=L.rASA[list(lidx)]
        rASA=R.rASA[list(ridx)]
        A.extend(lASA)
        A.extend(rASA)
        S=S+len(L.rASA)*len(R.rASA)
        Sr=Sr+np.nansum(L.rASA>0.05)*np.nansum(R.rASA>0.05)
    except:
        pdb.set_trace()
        continue
"""
#for (idx,f) in enumerate(fnames):
#    L=myPDB.loader(f)
#    cid=f[-16:-12]
#    lr=f[-11]!='l'
#    print np.nansum(L.ASA>0.05),len(L.ASA)
#    pdb.set_trace()
    
    rname=L.name.split('_')
    rname=rname[0]+'_r_'+rname[2]
    rname=pklpath+rname+'.pdb.pkl'
    
    R=myPDB.loader(rname)
    for c in R.stx[0]:
        rcid.append(c.id)
    cc=list(set(lcid)&set(rcid))
    if len(cc)>0:
        pdb.set_trace()
    
    #L.getStride()
    """
"""
    rASA.append(L.rASA)
    ASA.append(L.ASA)
    RD.append(L.RD[0])
    idxx=idx=np.where(~(np.isnan(L.RD[0])))
    nRD.append(L.RD[0]/L.RD[0][idx].max())
    CN.append(L.CN)
    UN.append(L.UN)
    DN.append(L.DN)
    PHI.append(L.Phi)
    PSI.append(L.Psi)    
    Z=np.nonzero(L.SS>0)
    X=np.empty(L.SS.shape[1]); X.fill(np.nan)
    X[Z[1]]=Z[0]
    SS.append(X)
ZrASA=np.concatenate(rASA)
idx=np.where(~(np.isnan(ZrASA)))
ZrASA=ZrASA[idx]
ZASA=np.concatenate(ASA)[idx]

ZRD=np.concatenate(RD)[idx]
ZnRD=np.concatenate(nRD)[idx]

ZCN=np.concatenate(CN)[idx]
ZUN=np.concatenate(UN)[idx]
ZDN=np.concatenate(DN)[idx]
ZPHI=np.concatenate(PHI)[idx]
ZPSI=np.concatenate(PSI)[idx]
ZSS=np.concatenate(SS)[idx]
#np.where(np.isnan(rASA))
cc=('r.','g.','b.','c.','m.','y.','k.','mx')
for k in range(6):
    idx=ZSS==k
    plt.plot(ZPHI[idx],ZPSI[idx],cc[k])
plt.show()
"""