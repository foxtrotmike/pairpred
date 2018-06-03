# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 19:50:07 2013
Code for displaying and analyzing the results from PAIRPred
@author: root
Usage:
    python analyzePredFile.py -i ./Flaviviridae/NS3_NS5_DBD4.InterPred.txt -o tmp.txt -n 10 -w 5 -l A -r A -t 2

"""
import numpy as np
#import matplotlib
#matplotlib.use('wx')
import matplotlib.pyplot as plt
from copy import deepcopy
import re
import sys
from optparse import OptionParser
# import pdb

def write2File(L,R,Scx,rZ,pwKfname,tx,auc=None,fname=None,lname='L',rname='R'):
    """
    Write the results to a file
    Format: (Example)
        AUC = 0.614399446175
        L_A R 258 68 R_A D 129 43 -0.307060015218 1
        L_A R 261 71 R_A D 129 43 -0.50449226754 1
        ...
        L_chainid <> res_name <> residue id in pdb file <> residue number in sequence or FASTA file <> residue number in L.R <>\
        R_chainid <> res_name <> residue id in pdb file <> residue number in sequence or FASTA file <> residue number in R.R <>\
        decision value <> label 
    AUC is only printed if evalROC was set to true
    labels are true labels if evalROC is true (otherwise all are -1)
    """
    lresi=dict (zip(L.resi.values(),L.resi.keys()))
    lR2Si=dict (zip(L.S2Ri,range(len(L.S2Ri))))
    rresi=dict (zip(R.resi.values(),R.resi.keys()))
    rR2Si=dict (zip(R.S2Ri,range(len(R.S2Ri))))
    lname=lname+':'
    rname=rname+':'
    #pdb.set_trace()
    if fname is None:
        fname='#'.join(Scx[0][0])+'.pairpred.txt'
    f = open(fname, 'w+')
    try:
        s='# Complex Name: '+str(Scx[0][0]);f.write(s+'\n')        
        s='# ligand file: '+L.ifname;f.write(s+'\n')
        s='# Receptor file: '+R.ifname;f.write(s+'\n')
        s='# Ligand Name: '+lname;f.write(s+'\n')       
        s='# Receptor Name: '+rname;f.write(s+'\n')       
        s='# Time taken (s): '+str(tx);f.write(s+'\n')
        s='# Number of examples: '+str(len(Scx));f.write(s+'\n')
        if auc is not None:           
            if type(auc)==type(tuple()):
                axx=''
                for a in auc:
                    axx+=str(a)+' '
                auc=axx
            s='AUC = '+str(auc)+'\n'
            f.write(s)
        for (k,(cid,(rl,rr),rlbl)) in enumerate(Scx):
            rla=rra='-'
            rlai=rrai=-1
            if rl in lR2Si:
                rla=L.seq[lR2Si[rl]]
                rlai=lR2Si[rl]
            if rr in rR2Si:
                rra=R.seq[rR2Si[rr]]
                rrai=rR2Si[rr]
            s=lname+lresi[rl][0]+' '+rla+' '+lresi[rl][1]+' '+str(rlai)+' '+str(rl)+' '+\
             rname+rresi[rr][0]+' '+rra+' '+rresi[rr][1]+' '+str(rrai)+' '+str(rr)+\
             ' '+str(rZ[k])+' '+str(rlbl)+'\n'
            f.write(s)
        
    except Exception as e:
        print "Error saving Results: ",e
        raise e
    finally:
        f.close()



def parseLine(lns,usePDBidx=True):
    if usePDBidx:
        xlid=2
        xrid=7
    else:
        xlid=4
        xrid=9
    lpid=lns[0].split(':')
    if len(lpid)<2:
        lcid=''
    else:
        lcid=lpid[1]    
    lid=int(re.sub("\D", "", lns[xlid])) #to remove het
    rpid=lns[5].split(':')
    if len(rpid)<2:
        rcid=''
    else:
        rcid=rpid[1]    
    rid=int(re.sub("\D", "", lns[xrid])) #re.sub("\D", "", "aas30dsa20")
    lrid=int(lns[4])
    rrid=int(lns[9])
    v=float(lns[10])
    l=int(lns[11])
    ls=lns[1]
    rs=lns[6]
    try:
        lpdbid=int(re.sub("\D", "", lns[2]))
        rpdbid=int(re.sub("\D", "", lns[7])) #[a-zA-Z_]
    except:
        pdb.set_trace()
    return (lcid,lid,rcid,rid,v,l,ls,rs,lrid,rrid,lpdbid,rpdbid)
    
def readFile(ifile,ilcid=None,ircid=None,usePDBidx=True,auconly=False):
    """
    ifile: Input result file
    ilcid: ligand chain id
    ircid: Receptor chain id
    usePDBidx: Whether to report the indices in the PDB file or in the 
        L.R and R.R
    auconly: When true, only the AUC is returned
    Both ilcid and ircid must be specified or both be set to None (equiv. to 
    leaving unspecified)
    Read result file (ifile) produced by PAIRPred and return (auc,Mv,Ml,lseq,rseq,lrV,rrV)
        auc: AUC score -- None if it does not exist in the file
        Mv: Matrix of scores with the rows corresponding to the ligand and 
            the columns corresponding to the receptor. 
        Ml: Matrix of labels with the rows corresponding to the ligand and 
            the columns corresponding to the receptor.             
        lseq: ligand sequence
        rseq: Receptor sequence
        lrV: Ligand scores (max)
        rrV: Receptor scores (max)
    """
    auc=None
    Lid=[]
    Rid=[]
    V=[]
    L=[]
    lseq={}
    rseq={}
    lrV={}
    rrV={}
    for ln in open(ifile,'r'):
        lns=ln.split()        
        if lns[0]=='#': #ignore comments
            continue
        elif lns[0]=='AUC':            
            auc=float(lns[2])
            if auconly:
                return auc
        else:
            (lcid,lid,rcid,rid,v,l,ls,rs,lrid,rrid,lpdbid,rpdbid)=parseLine(lns,usePDBidx)
            
            if lrid not in lrV:
                lrV[lrid]=(v,l)                
            else:
                lrV[lrid]=(np.max((lrV[lrid][0],v)),np.max((lrV[lrid][1],l)))
            if rrid not in rrV:
                rrV[rrid]=(v,l)                
            else:
                rrV[rrid]=(np.max((rrV[rrid][0],v)),np.max((rrV[rrid][1],l)))    
      
            if (ilcid is None and ircid is None) or (lcid==ilcid and rcid==ircid):
                Lid.append(lid)
                Rid.append(rid)
                #Lpdbid.append(lpdbid)
                #Rpdbid.append(rpdbid)
                V.append(v)
                L.append(l)
                lseq[lpdbid]=ls
                rseq[rpdbid]=rs
    #
    
    for ri in range(np.max(lseq.keys())):
        if ri not in lseq:
            lseq[ri]='-'
    for ci in range(np.max(rseq.keys())):
        if ci not in rseq:
            rseq[ci]='-'    
    
    assert len(np.unique(np.diff(np.sort(lseq.keys()))))==1 and (np.unique(np.diff(np.sort(lseq.keys()))))[0]==1
    assert len(np.unique(np.diff(np.sort(rseq.keys()))))==1 and (np.unique(np.diff(np.sort(rseq.keys()))))[0]==1
    Mv=np.zeros((np.max(Lid)+1,np.max(Rid)+1))
    Mv.fill(np.nan)
    Ml=np.zeros(Mv.shape)
    Ml.fill(np.nan)
    #pdb.set_trace()
    #lseq=lseq.values();lseq=''.join(lseq)
    #rseq=rseq.values();rseq=''.join(rseq)     
    for i in range(len(Lid)):
        if Lid[i]>=0 and Rid[i]>=0:
            Mv[Lid[i],Rid[i]]=np.nanmax((Mv[Lid[i],Rid[i]],V[i]))
            Ml[Lid[i],Rid[i]]=np.nanmax((Ml[Lid[i],Rid[i]],L[i]))
    return (auc,Mv,Ml,lseq,rseq,lrV,rrV)
def M2LR(Mvc):
    """
    convert complex wise scores to individual protein scores
    """
    Lv0=np.nanmax(Mvc,axis=1)
    Rv0=np.nanmax(Mvc,axis=0)
    return Lv0,Rv0
def rV2V(lrV):
    """
    Return the prediction scores from dictionary object lrV
    """
    return np.array([a for (a,b) in lrV.values()])
def plotPred(Mo,lrV,rrV,thr=2.5,lname='L',lseq='',rname='R',rseq='',met='PAIRPred'):
    """
    Plots the results from PAIRPred. It displays the matrix of scores as an image
    It also generates a imahe in which only the top 'thr' % of scores are displayed.
    It generates the plots for scores for the ligand and receptor proteins as well.
        Mo: score matrix
        thr: select top thr % only
        lname: name of the ligand
        lseq: sequence of the ligans
        rname: ma,e pf tje receptor
        rseq: Seqiemce of the receptor
        met: Name of the method
    """
    xstep=15
    ystep=15
    met=met+': '
    cmap =plt.cm.jet;# plt.get_cmap()
    cmap.set_bad(color = 'k', alpha = 1.)
    Mf=Mo.flatten()
    Mf=Mf[~np.isnan(Mf)]
    plt.figure()
    (c,b,p)=plt.hist(Mf,bins=1000,cumulative=True,normed=True);#plt.xlabel('value');plt.ylabel('Normalized Frequency');plt.title(met+'distribution')
    plt.close()
    thrv=b[np.where(c<(1-thr/100.0))[0][-1]]
    #xtn=range(0,Mo.shape[1],50);xt=[Lpid[x] for x in xtn]
    plt.figure();plt.imshow(Mo,cmap=cmap);plt.colorbar();plt.xlabel(rname);plt.ylabel(lname);plt.title(met+'All Predictions');plt.xticks(range(0,Mo.shape[1],xstep),rotation='vertical');plt.yticks(range(0,Mo.shape[0],ystep));
    M1=deepcopy(Mo)
    M1[Mo<thrv]=np.nan
    if type(rseq)==type(dict()):
        ridx=[pi for pi in rseq.keys() if rseq[pi]!='-']
    else:
        ridx=range(0,Mvc.shape[1])
    if type(lseq)==type(dict()):        
        lidx=[pi for pi in lseq.keys() if lseq[pi]!='-']
    else:
        lidx=range(0,Mvc.shape[0])
    #pdb.set_trace()
    Lv=rV2V(lrV)
    Rv=rV2V(rrV)
    plt.figure();plt.imshow((M1),cmap=cmap);plt.colorbar();plt.xlabel(rname);plt.ylabel(lname);plt.title(met+'Top '+str(thr)+'% Predictions');plt.xticks(range(0,Mvc.shape[1],xstep),rotation='vertical');plt.yticks(range(0,Mvc.shape[0],ystep));
    plt.figure();plt.stem(ridx,Rv,'-.');plt.title(met+rname+' Prediction');plt.xticks(range(np.min(ridx),np.max(ridx),xstep),rotation='vertical');plt.grid()
    plt.figure();plt.stem(lidx,Lv,'-.');plt.title(met+lname+' Prediction');plt.xticks(range(np.min(lidx),np.max(lidx),ystep),rotation='vertical');plt.grid()

def sortScores(Mo):
    M1=deepcopy(Mo)
    Mor=M1.ravel()
    nidx=np.isnan(Mor)
    Mor[nidx]=-np.inf
    Mori=np.argsort(Mor)[::-1]
    Mori=[i for i in Mori if Mor[i]!=-np.inf]
    (r,c)=np.unravel_index(Mori,Mo.shape)
    v=Mor[Mori]
    
    return (r,c,v)
    
def writeTop2File(Mo,fname,lseq='',rseq='',HWS=5,top=np.inf):
    """
    Given the score matrix, it writes the top 'top' pairs into a file
    178<\t>	247<\t>1.43149930273<\t>EPDYEVDEDIF<\t>	V<\t>	RFTMRHKKATY<\t>H
    ligand_resid,receptor_resid,score,HWS sized sequnce window around the residue, residue...
    """
    
    def getslseq(lseq,HWS,ri):
        rle=ri-HWS
        rre=ri+HWS+1
        rlpad=rrpad=''
        if rle<1:                
            rlpad='*'*(0-rle+1)
           #pdb.set_trace()
            rle=1            
        if rre>len(lseq)-1:
            rrpad='*'*(rre-len(lseq)+1)
            #pdb.set_trace()
            rre=len(lseq)-1            
        slseq=rlpad+lseq[rle:rre]+rrpad
        return slseq    
        
    (r,c,v)=sortScores(Mo)
    f = open(fname, 'w+')
    s='#lig_id\trec_id\tscore\tlig_res\tlig_win\trec_res\trec_win\n'
    f.write(s)
    for i in range(int(np.min((len(v),top)))):
        s=str(r[i]) + '\t'+str(c[i])+'\t'+str(v[i])
        if len(lseq):
            slseq=getslseq(lseq,HWS,r[i])
            s=s+'\t'+lseq[r[i]]+'\t'+slseq
        else:
            s=s+'\t'+'-'
        s=s+'\t'
        if len(rseq):
            srseq=getslseq(rseq,HWS,c[i])
            s=s+'\t'+rseq[c[i]]+'\t'+srseq
        else:
            s=s+'\t'+'-'      
        
        f.write(s+'\n')
    f.close()
    #return (r,c,v)
if __name__=="__main__":
    #ifile,ilcid,ircid,topfname,top,HWS,plotthr
    parser = OptionParser()
    parser.add_option("-i", "--ifile", dest="ifile",
                      help="Input Prediction file", metavar="file path")
    parser.add_option("-l", "--lchain", dest="ilcid",
                      help="Ligand chain (default: All)", metavar="character")
    parser.add_option("-L", "--lname", dest="lname",
                      help="Ligand name (default: L)", metavar="character")         
    parser.add_option("-R", "--rname", dest="rname",
                      help="Receptor name (default: R)", metavar="character")                          
    parser.add_option("-r", "--rchain", dest="ircid",
                      help="Receptor chain (default: All)", metavar="character")
    parser.add_option("-o", "--ofile", dest="ofile",
                      help="Output File", metavar="file path")
    parser.add_option("-n", "--num", dest="num",
                      help="number of top predictions in output file (default: 100)", metavar="Number")    
    parser.add_option("-w", "--win", dest="win",
                      help="sequence half-window size (default: 5)", metavar="Number")             
    parser.add_option("-t", "--thr", dest="thr",
                      help="top % threshold in plotting (default: 1 %)", metavar="Percentage number")   
    parser.add_option("-p", "--pp", dest="post",
                      help="Whether to perform postprocessing or not. When empty or unsepecified, no post-processing otherwise it should point to the path of the cid_l_u and cid_r_u pdb files", metavar="string")                        
    (options, args) = parser.parse_args()
    
    if options.ifile is None or options.ofile is None:
        print 'Must specify both input and output file names. Exiting.'
        print 'Type python analyzePredFile.py -h for help'
        sys.exit(1)
    if options.post is None or len(options.post)==0:
        post=False
    else:
        post=True
        from postProcess import *
        pp_cid=getFileParts(getFileParts(options.ifile)[1])[1]
        pp_pdbpath=os.path.join(options.post,'')
        pp_ppath=os.path.join(getFileParts(options.ifile)[0],'')
    
    if options.lname is None:
        options.lname='L'
    if options.rname is None:
        options.rname='R'        
    if options.win is None:
        options.win=5
    else:
        options.win=int(options.win)
    if options.thr is None:
        options.thr=1.0
    else:
        options.thr=float(options.thr)
    if options.num is None:
        options.num=100
    else:
        options.num=int(options.num)    
    print options
    if not post:
        (auc,Mvc,Mlc,lseq,rseq,lrV,rrV)=readFile(options.ifile,ilcid=options.ilcid,ircid=options.ircid)
    else:
        auc,Mvc0,Mvc,Mlc,lseq,rseq,lrV0,lrV,rrV0,rrV=postProcessAvg(cid=pp_cid,pdbpklpath=pp_pdbpath,pppath=pp_ppath)
        #(pauc,Mv,Ml,lseq,rseq,lauc,rauc)
        from analyzeDists import *
        Z=computeDistMeansForComplex(cid=pp_cid,N=100,pdbpklpath=pp_pdbpath,pppath=(auc,Mvc,Mlc,lseq,rseq,lrV,rrV))
        print Z
    if auc is not None:
        print 'AUC =', auc
    
    #writeTop2File(Mvc,fname=options.ofile,lseq=lseq,rseq=rseq,HWS=options.win,top=options.num)
    plotPred(Mvc,lrV,rrV,lname=options.lname,rname=options.rname,met='PAIRPred',thr=options.thr,lseq=lseq,rseq=rseq)
    plt.show()