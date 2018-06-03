# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 21:55:38 2012
@author: Afsar
This module constructs positive and negative training examples from a data set of pdb files 
"""
import glob
import os
import cPickle
from getPair import *
import itertools
import random
import numpy as np

def getDist2NearestEx(P,Npure,dLL,dRR):
    pd=np.zeros(len(P))
    for i,p in enumerate(P):
        pd[i]=np.inf
        for n in Npure:
            d=dLL[p[0],n[0]]+dRR[p[1],n[1]]
            if d<pd[i]:
                pd[i]=d
    return list(pd)


def getPosex(idir,cname,pdthr=6.0,mdthr=8.0):
    luname=os.path.join(idir,cname[0]+'.pdb')
    runame=os.path.join(idir,cname[1]+'.pdb')
    lbname=os.path.join(idir,cname+'_l_b.pdb')
    rbname=os.path.join(idir,cname+'_r_b.pdb')
    C=myPDBComplex([lbname,rbname],mdthr )
    (D,LenU)=C.findNforU([luname,runame])
    #pdb.set_trace()
    D=D[0][1]
    lenL=LenU[0]
    lenR=LenU[1]
    dict_gdf=dict()
    P=[]
    Pl=set([])
    Pr=set([])
    for (i,j,d) in D:
        if ~ (np.isnan(i) or np.isnan(j)):                      
            dict_gdf[(i,j)]=d
            if d<=pdthr: # positive examples
                Pl.add(i)
                Pr.add(j)
                P.append((i,j))
            else:
                pass #residue pairs whose distance is > pdthr  but less than dthr
                #keep a list of those and their distances and once
                #negative examples have been selected, see if the
                #residue pair occurs in the list or not. If it does
                # then calculate the distance, otherwise no    
    return (Pl,Pr,P,lenL,lenR)
class getExamplesDBD:
    """
    An object of this class will contain two dictionaries Pex amd Nex with 
    the same keys (complex_ids)
    Pex[key][0] is a list of tuples (L.R_idx,R.R_idx) of indicies from L and R
        proteins for complex id key that lie within a distance of 'dthr' of one 
        another. 
    Pex[key][1] is the number of residues in L
    Pex[key][2] is the number of residues in R
    Nex[key] is a list of tuples of indices from L and R that constitute negative
        examples
    dir: is the directory in which the pdb files were located from which the
        object has been constructed
    pdthr: distance threshold for positive examples
    mdthr: distances are saved in dict_gdf if they are less than this threshold
    dict_gdf[cname][key]: is a dictionary object with key = lres,rres) and the distance between them as the value, not saved for every example
    dpd[cname][key]: is a dictionary object with key=(lres,rres) that saves distance of the positive example key  to its nearest pure negative example
    dnd[cname][key]: is a dictionary object with key=(lres,rres) that saves distance of the negative example key  to its nearest posiive example
    Note: If there is an error in processing any file in a complex, the complex
        is not included in the object
    """
    def __init__(self,idir,pfile):
        """
        idir: path to the directory in which the files are located
        pdthr: distance threshold (default: 6.0 Angstroms) for positive ex.
        mdthr: distances are saved for upto mdthr in dict_gdf
        """        
        self.P=myPickle.load(pfile)
        self.Pex=dict()
        self.Nex=dict()
        self.Plr=dict()
        
        for cid in self.P:      
            P=self.P[cid]
            for (lc,rc) in P:
                cname=(cid+'_l_u~'+lc.strip(),cid+'_r_u~'+rc.strip())
                print "############################ Processing:",cname
                #try:                
                (_,lR,_,_,_)=readPDB(os.path.join(idir,cname[0].strip()+'.pdb'))
                (_,rR,_,_,_)=readPDB(os.path.join(idir,cname[1].strip()+'.pdb'))
                lenL=len(lR)
                lenR=len(rR)
                Pc=np.asarray(P[(lc,rc)],int)
                Pl=list(np.unique(Pc[:,0]))
                Pr=list(np.unique(Pc[:,1]))                
                Pc=[tuple(x) for x in Pc]
                
                    #(Pl,Pr,Pc,lenL,lenR)=getPosex(idir,cname)                        
#                except:                
#                    print "!!!!!!!!!!!!!!!!!!!!!!!! Error processing: "+cname
#                    continue
                
                self.Nex[cname]=[]
                self.Pex[cname]=(Pc,lenL,lenR)
                self.Plr[cname]=[Pl,Pr]
    def __repr__(self):
        s='getExamplesDBD Instance'   
        s+='\n\t Number of complexes: '+str(len(self.Pex))
        s+='\n\t Number of positive examples: '+str(np.sum([len(self.Pex[k][0]) for k in self.Pex]))
        s+='\n\t Number of Negative examples: '+str(np.sum([len(self.Nex[k]) for k in self.Nex]))       
        return s
    def save(self,ofname=None):
        """
        Save the object
        """
        if ofname is None:
            ofname='E_'+str(self.pdthr)+'.lbl.pkl'
        print "File Saved: "+ofname
        output = open(ofname, 'wb')
        cPickle.dump(self, output)
        output.close()
        
    @classmethod   
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """
        return cPickle.load(open(pklfile, "rb" ) )
    def getPureNegEx(self,cname):
        """
        Return all negative examples in which both residues do not interact with any residue in the other protein
        """        
        return [(i,j) for (i,j) in self.getNegEx(cname) if i not in self.Plr[cname][0] and j not in self.Plr[cname][1]]
    """
    def attachD2N(self):
        
        #Attach distances to nearest example of opposite sign for both positive (dpd) and negative (dnd) examples
        
        self.dpd={}
        self.dnd={}
        for cname in E.Pex.keys():
            print "Processing :",cname
            (_,lres,_,_,_)=readPDB(os.path.join(self.dir,cname+'_l_u.pdb'))            
            (_,rres,_,_,_)=readPDB(os.path.join(self.dir,cname+'_r_u.pdb'))                 
            dLL=getDistMat(getCoords(lres))
            dRR=getDistMat(getCoords(rres))
            Npure=self.getPureNegEx(cname)
            N=self.Nex[cname]
            P=self.Pex[cname][0]
            pd=getDist2NearestEx(P,Npure,dLL,dRR)
            nd=getDist2NearestEx(N,P,dLL,dRR)
            self.dpd[cname]=dict(zip(P,pd))
            self.dnd[cname]=dict(zip(N,nd))
    """            
    def getNegExSingle(self,cname,p=1,mode='all'):
        Nl=list(set(range(self.Pex[cname][1])).difference(self.Plr[cname][0]))
        Nr=list(set(range(self.Pex[cname][2])).difference(self.Plr[cname][1]))
        if mode=='+': #same number as the number of positive examples
            Nl=random.sample(Nl,np.min((len(Nl),p*len(self.Plr[cname][0]))))
            Nr=random.sample(Nr,np.min((len(Nr),p*len(self.Plr[cname][1]))))   
        return (Nl,Nr)
        
    def getNegEx(self,cname):
        """
        Returns all negative examples for a given complex 'cname'
        LxR - P
        """
        try:
            L=[range(self.Pex[cname][1]),range(self.Pex[cname][2])]
            X=set(itertools.product(*L))
            return X.difference(set(self.Pex[cname][0]))
        except:
            return None
            
    def selectNegEx(self,cname=None,p=100.0,M=1000):    
        """
        Samples negative examples
        cname: list of complex ids for which the sampling is to be done
            if not specified, then the function automatically performs the
            sampling for all complexes in the object
        M: selection of sampling technique
            if M is not specified, then M=1000 and the minimum of M and 
                the number of negative examples in a complex is taken
            if M is set to None, then all negative examples are included
            If M is '+-' then the negative examples is equal to the minimum of 
                the p (default 100) % of the number of positive examples and 
                the number of all possible negative examples in a complex and the
                ngative examples are selected as follows:
                    At most 20% of negative examples come from pairs of residues each 
                    of which interacts with some other residue but not with each other  
                    At most 30% of negative examples come from pairs of residues exactly one of
                    which is interacting with some other residue
                    The remaining examples are "true negative" examples, i.e., neither of the
                    residues interacts with any other residue
                    
            if M is '+' then the number of negative examples is taken to 
                be equal to the minimum of the p (default 100) % of the number of positive examples and 
                the number of all possible negative examples
            Otherwise, M is expected to be a number and the number of negative
                examples is taken to be the minimum of p% of the total number
                of negative examples and M
            Negative examples are then randomly sampled and stored in Nex
        """
        if cname is None:
            cname=self.Pex.keys()
        if M is None:
            M=np.inf            
        for c in cname:            
            if M=='+-':
                Pl,Pr=self.Plr[c] #Positive examples in L and R
                Nl=range(self.Pex[c][1]) #All examples in L
                Nr=range(self.Pex[c][2]) #All examples in R                
                XPlPr=set(itertools.product(Pl,Pr)).difference(set(self.Pex[c][0])) #Negative examples from (+,+) pairs
                XPlNr=(set(itertools.product(Pl,Nr)).difference(set(self.Pex[c][0]))).difference(XPlPr) #Negative examples from (+,-) pairs
                XNlPr=(set(itertools.product(Nl,Pr)).difference(set(self.Pex[c][0]))).difference(XPlPr) #Negative examples from (-,+) pairs
                XPN=XPlNr.union(XNlPr) #negative examples involving exactly 1 positive and 1 negative, i.e., (+,-) or (-,+)
                XNlNr=set(itertools.product(Nl,Nr)).difference(set(self.Pex[c][0])) #Negative examples from (-,-)
                XNlNr=XNlNr.difference(XPlPr)
                XNlNr=XNlNr.difference(XPN)
                npx=int(p*len(self.Pex[c][0])/100.0)
                #(+,+): 20%, (+,-) or (-,+):30%, (-,-): 50%
                npp=int(min(0.15*npx,len(XPlPr)))
                N=set(random.sample(XPlPr,npp))
                npn=int(min(0.35*npx,len(XPN)))
                N=N.union(set(random.sample(XPN,npn)))
                nn=npx-len(N)
                #pdb.set_trace()
                N=N.union(set(random.sample(XNlNr,nn))) # was originally mistakingly XPN before 22APR2013
                self.Nex[c]=list(N)                        
            else:
                N=self.getNegEx(c)
                if M=='+':
                    n=int(min(p*len(self.Pex[c][0])/100.0,len(N)))
                else:
                    n=int(min(np.floor(p*len(N)/100.0),M))
                self.Nex[c]=random.sample(N,n)
if __name__=="__main__":    
    pdbpath= '../../DBD4CSPKL/PDB_' #'/s/chopin/b/grad/minhas/Desktop/DBD4N/DBD4/'#'./DBD4N/DBD4/'#
    E=getExamplesDBD(pdbpath,pfile='../../DBD4CSPKL/DBD4_.pos.mkl')
    
    print "-----------Choosing Negative Examples-----------"
    E.selectNegEx(p=125.0,M='+-')
    print "---Computing distances to nearest examples of opposite sign----"
    #E.attachD2N()
    E.save('test_.lbl.pkl')
            
            