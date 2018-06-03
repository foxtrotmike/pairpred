# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:03:11 2013

@author: root
"""

def attachD2N(self,pdbpklpath):
    """
    Attach distances to nearest positive or negative examples
    """
    self.dpd={}
    self.dnd={}
    for cname in E.Pex.keys():
        print "Processing :",cname
        L=myPDB.loader(os.path.join(pdbpklpath,cname+'_l_u.pdb.pkl'))
        R=myPDB.loader(os.path.join(pdbpklpath,cname+'_r_u.pdb.pkl'))
        dLL=getDistMat(L.Coords)
        dRR=getDistMat(R.Coords)
        Npure=getPureNegEx(E,cname)
        N=E.Nex[cname]
        P=E.Pex[cname][0]
        pd=getDist2NearestEx(P,Npure)
        nd=getDist2NearestEx(N,P)
        self.dpd[cname]=dict(zip(P,pd))
        self.dnd[cname]=dict(zip(N,nd))

pdbpklpath='./DBD4N/PDBPKL4/'
