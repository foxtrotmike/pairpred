# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 13:41:30 2013
For analysising the results from NS3-NS5 predictions
@author: root
"""
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
# import pdb
from testSingleComplex_par import readFile,plotPRED

    
if __name__=="__main__":
    ifile='../Flaviviridae/1L9K.InterPRed.txt'   
    (auc,Mo,lseq,rseq)=readFile(ifile)
    plotPRED(Mo,lname='1L9K',rname='1L9K',met='PAIRPred-Struct',thr=1.0)
    plt.show()
    """
    Ms=np.zeros((618,900))
    Ms.fill(np.nan)
    
    sfile='shandarpred.txt'
    for ln in open(sfile):
        lns=ln.split()
        Ms[int(lns[0][1:]),int(lns[1][1:])]=float(lns[2])
        #pdb.set_trace()
    Ms[np.isnan(Ms)]=np.nanmin(Ms)
    
    
    #plotPRED(Ms,met='PPIPP')
    plt.show()
    """