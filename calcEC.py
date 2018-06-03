
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:02:07 2013
Wrapper for calcEC.m to compute the mutual and direct information from a MSA
@author: root
"""

import os
import tempfile
# import pdb
import sys
mlabpath='/s/parsons/l/sys/matlab64/bin/'
if mlabpath not in  os.environ["PATH"]:
    print 'Adding '+mlabpath+' to system path'
    os.environ["PATH"]+=os.pathsep+mlabpath
    
def parseECFile(fpath):
    """
    Given a file created by calcEC.m, returns dictionaries for MI and DI in 
    which the key is the tuple (i,j) indicating the indices of the pairse of
    locations in the sequence, i>j
    """
    M={}
    D={}
    for f in open(fpath,'r'):
        f=f.split()
        r=int(f[0])
        c=int(f[2])
        mi=float(f[4])
        di=float(f[5])
        M[(r-1,c-1)]=mi
        D[(r-1,c-1)]=di
    return M,D
      
def runEC(msa,cid,ofile=None,mlab='matlab',pgm='calcEC'):
    """
    msa: single msa file in fasta format, comma seperated string list of such files
        or list of msa files
    cid: sequence id of interest in the first msa file
    ofile (optional): When specified the output of calcEC.m is saved here
    mlab: matlab path
    pgm: path to calcEC.m
    Return (same as for parseECfile)
    """
    msal=msa  
    if type(msa)==type([]):
        msal=','.join(msa)
    msa='{'+msal+'}'
    if ofile is None:
        out_file = tempfile.NamedTemporaryFile(suffix='.ec')
        out_file.close()
        ofname=out_file.name
    else:
        ofname=ofile
    cmdstr=mlab+' -nosplash -nodesktop -nodisplay -nojvm -wait -r \"'+pgm+' '+ msa +' ' +cid+ ' '+ofname+'; exit(0)\"'
    #cmdstr=mlab+' -nosplash -nodesktop -nodisplay -nojvm -r \''+pgm+' '+ msa +' ' +cid+ ' '+ofname+'; exit(0)\''
    #pdb.set_trace()
    if os.system(cmdstr):
        print "Error running",cmdstr
        R=[]
    else:        
        R=parseECFile(ofname)
        if ofile is None:
            os.remove(ofname)
    return R
    
if __name__=="__main__":
    usage="python calcEC.py msafile cid ofile"
    if len(sys.argv)>1:
        msa=sys.argv[1]
        cid=sys.argv[2]       
        ofile=sys.argv[3]
        runEC(msa,cid,ofile)
    else:
        print "The usage is:", usage