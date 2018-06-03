# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:17:59 2013

@author: root
"""

import Bio.pairwise2
from myPDB import *
from Bio.pairwise2 import format_alignment 

bdir='../../DBD4N/PDBPKL4/'
cid='2C0L'
Lu=myPDB.loader(bdir+cid+'_l_u.pdb.pkl')
Lb=myPDB.loader(bdir+cid+'_l_b.pdb.pkl')

Ru=myPDB.loader(bdir+cid+'_r_u.pdb.pkl')
Rb=myPDB.loader(bdir+cid+'_r_b.pdb.pkl')

us=Lu.seq
bs=Lb.seq

aln=Bio.pairwise2.align.globalxs(us,bs,-1,-0.1)[0]

print format_alignment(*aln) 