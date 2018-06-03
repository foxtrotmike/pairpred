# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 20:26:14 2013
Replication of "Relative Solvent Accessible Surface Area Predicts Protein 
Conformational Changes upon Binding" for DBD 4.
Saves rmsd_atomic.mkl file which contains a dictionary:
    cid_l/r: rmsd, unbound_asa, bound_asa, unbound_mol_weight,bound_mol_weight, atomic_rmsd_unbound, atomic_rmsd_bound
@author: root
"""

from pymol import cmd, CmdException

import Bio.PDB
import numpy as np
from BISEPutils import readPDB,getFileParts
import Bio.pairwise2
from getPAIR import mapU2B
from myPDB import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Polypeptide import three_to_one
# import pdb
import myPickle
import os.path
import bottleneck as BN
from dbd4 import parseCSVData
import time
#import random
def getBvalues(R):
    """
    Save the b values for each residue        
    """
    B=[[] for _ in range(len(R))]
    for (idx,r) in enumerate(R):
        #pdb.set_trace()
        B[idx]=[ak.get_bfactor() for ak in r ]
    return B
    
def calcRMSD_pymol(uf,bf):
    """
    Given two pdb files of the same protein, this function calculates the
    rmsd, asa for each and the molecular weight of each using Pymol
    """        
    # Call the function below before using any PyMOL modules.
    #time.sleep(random.random())
    cmd.set("dot_solvent", 1)
    cmd.load(uf)
    cmd.load(bf)
    #cmd.h_add()
    #cmd.remove('het')
    _,un,_=getFileParts(uf)
    _,bn,_=getFileParts(bf)
    asa_u=cmd.get_area(un)
    asa_b=cmd.get_area(bn)
    
    umass=cmd.get_model(un).get_mass()
    bmass=cmd.get_model(bn).get_mass()
    #rms=cmd.super(un,bn,transform=1)[0]    
    #time.sleep(random.random())
    bv0=[]
    cmd.iterate( 'all', 'bv0.append(b)', space=locals())
    cmd.do('run colorbyrmsd.py; colorbyrmsd \''+un+'\',\''+bn+'\',guide = 0,doAlign=1, doPretty=1')
    while True:        # synchronization
        bv1=[]
        cmd.iterate( 'all', 'bv1.append(b)', space=locals())                
        if bv0!=bv1:
            time.sleep(0.1)
            break
        
    out_file = tempfile.NamedTemporaryFile(suffix='.pdb'); out_file.close();tmp_pdb=out_file.name
    updb=tmp_pdb+'u';bpdb=tmp_pdb+'b'    
    cmd.save(updb,un);cmd.save(bpdb,bn)    
    (_,uR,_,_,_)=readPDB(updb); urmsd=getBvalues(uR);os.remove(updb)
    (_,bR,_,_,_)=readPDB(bpdb); brmsd=getBvalues(bR);os.remove(bpdb)
    rms=np.sqrt(np.mean(np.array([v for V in urmsd for v in V if v>=0])**2))
    #(_,urx,_,_,_)=readPDB(uf); ux=getBvalues(urx);
#    if np.abs(rms-rmsd)>0.1:
#        print "RMSD =",rms,rmsd
#        pdb.set_trace()
    
    
    cmd.reinitialize()
    pdb.set_trace()
    return rms,asa_u,asa_b,umass,bmass, urmsd,brmsd
def calcRMSD(uf,bf):
    """
    RMSD using Biopython (not very good)
    """
    # ORiginal code derived from the ref below
    # Peter Cock, http://www.warwick.ac.uk/go/peter_cock/python/protein_superposition/
    if type(uf)==type(''):        
        (ustx,uR,_,useq,uS2Ri)=readPDB(uf)
    else:
        ustx,uR,useq,uS2Ri=uf.stx,uf.R,uf.seq,uf.S2Ri
    if type(bf)==type(''):             
        (bstx,bR,_,bseq,bS2Ri)=readPDB(bf)
    else:
        bstx,bR,bseq,bS2Ri=bf.stx,bf.R,bf.seq,bf.S2Ri    
    (u2b,b2u)=mapU2B(useq,uS2Ri,len(uR),bseq,bS2Ri,len(bR))
    """
    aln=Bio.pairwise2.align.globalxs(useq,bseq,-1,-0.1)[0]
    from Bio.Alphabet import IUPAC, Gapped
    z=Bio.Align.Generic.Alignment(Gapped(IUPAC.protein, "-"))
    z.add_sequence('1',aln[0])
    z.add_sequence('2',aln[1])
    import Bio.PDB.StructureAlignment as sa
    zz=sa(z,ustx[0],bstx[0])
    for (r1,r2) in zz.get_iterator():
        print (r1,r2)
    pdb.set_trace()
    """
    uA=[]
    bA=[]
    Uidx=[]
    Bidx=[]
    for uidx,bidx in enumerate(u2b):
        if ~np.isnan(bidx) and 'CA' in uR[uidx] and 'CA' in bR[bidx]:
            uA.append(uR[uidx]['CA'])
            bA.append(bR[bidx]['CA'])
            Uidx.append(uidx)
            Bidx.append(bidx)
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(uA, bA)
    super_imposer.apply(bstx[0].get_atoms())
    return super_imposer.rms,Uidx,Bidx
def getProps(f):
    """
    Code for getting the molecular weight and other properties using Biopython
    """
    L=myPDB.loader(f)
    aseq = ProteinAnalysis(L.seq)
    return aseq.molecular_weight(),np.max(aseq.flexibility()),np.sum(L.ASA)
def batchExtract(pkldir,bdir,ofname):
    """
    Running the information required for all files
    """
    import glob
    
    flist=glob.glob(pkldir+'*.pdb.pkl')
    TT=len(flist)+0.0
    if os.path.isfile(ofname) is False:
        fdict={}
    else:
        fdict= myPickle.load(ofname)  
    for cnt,f in enumerate(flist): 
        print '% Done =',cnt/TT
        (_,k,_)=getFileParts(getFileParts(f)[1])
        #pdb.set_trace()
        k=k[:-2]
        if k not in fdict:
            print "Processing",f
            try:
                U=myPDB.loader(pkldir+k+'_u.pdb.pkl')
                B=myPDB.loader(pkldir+k+'_b.pdb.pkl')
            except:
                continue      
            pdb.set_trace()
            #rmsd,Uidx,Bidx=calcRMSD(U,B)
            try:
                rpymol=calcRMSD_pymol(bdir+k+'_u.pdb',bdir+k+'_b.pdb')
            except:
                print "Error processing", k
                cmd.reinitialize()
                time.sleep(0.1)
                continue                
                
            #pdb.set_trace()
            #useq=''.join([three_to_one(U.R[i].get_resname()) for i in Uidx])
            #bseq=''.join([three_to_one(B.R[i].get_resname()) for i in Bidx])
            #a_useq=ProteinAnalysis(U.seq)
            #a_bseq=ProteinAnalysis(B.seq)
            #asa_u=np.sum([U.ASA[i] for i in Uidx])
            #asa_b=np.sum([B.ASA[i] for i in Bidx])
            fdict[k]=rpymol#+(BN.nanmean(U.B),BN.nanmean(B.B),BN.nanmedian(U.B),BN.nanmedian(B.B),BN.nanmax(U.B),BN.nanmax(B.B))
            #pdb.set_trace()
            myPickle.dump(ofname,fdict)
            print k,rpymol[0]  
        else:
            print "Already found",f            
    return fdict            
       

if __name__=="__main__":
    
#    cid='1Y64_r'
#    ddir='../DBD4N/DBD4/'
#    bf=ddir+cid+'_b.pdb'
#    uf=ddir+cid+'_u.pdb'
#    print "RMSD =",calcRMSD(uf,bf)
    ofname='rmsd_atomic.mkl'  
    loadDirect=True
    ddir='../DBD4N/DBD4/'
    pkldir='../DBD4N/PDBPKL4/'
    if not loadDirect or not os.path.isfile(ofname):
        import __main__
        __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI #
        import pymol
        pymol.finish_launching()
        from pymol import cmd
        fdict=batchExtract(pkldir,ddir,ofname)         
    else:
        fdict= myPickle.load(ofname)  
    #cid_l/r: rmsd, (rmsd_pymol, unbound_asa_pymol, bound_asa_pymol), unbound_mol_weight, bound_mol_weight,unbound_asa, bound_asa     
   
    Va=[]
    for v in fdict.values(): 
        Va.append(v[:-2])

    
    V=np.array(Va )   #rmsd, unbound_asa, bound_asa, unbound_mol_weight,bound_mol_weight
    
    #x=V[:,2]/(4.84*(V[:,4]**0.76)); y=V[:,1];  #x=x[nidx];y=y[nidx]; 
    x=V[:,1]/(4.84*(V[:,3]**0.76)); y=V[:,0];
    #x=V[:,3]/(0.346*V[:,5]+2.5e+03)
    #x=V[:,5]/(4.84*(V[:,2]**0.76)); y=V[:,1]; 
    nidx=(y>1e-3); 
    #nidx=nidx*(np.abs(V[:,4]-V[:,5])<1e3)
    x0=x[nidx];y0=y[nidx]; 
    
    from scipy.stats import spearmanr,pearsonr
    print "r value for relative ASA and over all rmsd at the protein level:",spearmanr(x0,y0),pearsonr(x0,np.log2(y0))
    
    dbd4=parseCSVData('..\Complete Data\DBD4_data.csv')
    #cid: cname, categ, lname, rname, rmsd, bsa
    dx={}
    for idx,key in enumerate(fdict):
        dx_key=key[:4]
        px=True and nidx[idx]
        if dx_key in dx:
            vx=(dx[dx_key][0]+x[idx])/2.0
            vy=(dx[dx_key][2]+y[idx])/2.0
#            if y[idx]>dx[dx_key][2]:
#                
#                px=True
#            else:
#                px=False
        else:
            vx=x[idx]
            vy=y[idx]
        if px:
            try:
                dx[dx_key]=vx,dbd4[dx_key][1],vy,dbd4[dx_key][-2],dbd4[dx_key][-1]
            except KeyError:
                continue
    (v,l,ormsd,irmsd,bsa)=zip(*dx.values())
    #from PyML.evaluators import roc
    #(_,_,auc)=roc.roc(list(v),list(l))
    print "r value for relative ASA and over all rmsd at the complex level:",spearmanr(v,ormsd)
    #print auc