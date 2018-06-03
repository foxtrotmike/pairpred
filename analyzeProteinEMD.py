# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 18:00:20 2013
Multidimensional Empirical Mode decomposition analysis of Protein properties
@author: Afsar
"""

from myPDB import *
from MEMD import *
from scipy.stats import nanmean
from getExamplesDBD import *
bdir='../DBD4N/PDBPKL4/'
cid='1A2K'
rl='r'
efile=bdir+'E_125PN_15_35_50.lbl.pkl'
E=getExamplesDBD.loader(efile)
L=myPDB.loader(bdir+cid+'_'+rl+'_u.pdb.pkl')
x0=np.vstack((L.ASA,L.RD[0],L.B))
x=x0
x=(x.T-nanmean(x,axis=1).T).T
x[np.isnan(x)]=0
x=x.T/np.std(x,axis=1)
M=MEMD(x);
B=np.unique([e[rl=='r'] for e in E.Pex[cid][0]])
print B
plt.figure();
for d in range(x.shape[1]):
    plt.subplot(1,x.shape[1],d+1)
    M.showHHSpectrum(d);
    ax=plt.gca()
    for b in B:
        ax.annotate('', xy=(b, 200), xytext=(b, 210),
                    arrowprops=dict(facecolor='black', shrink=0.0,width=0.5),
                    )
plt.show()