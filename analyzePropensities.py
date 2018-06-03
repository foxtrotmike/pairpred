# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 00:04:46 2013
Plot the propensities of different residues to bind to other residues
@author: root
"""
import os.path
from myPDB import *
from getExamplesDBD import *
from Bio.PDB.Polypeptide import three_to_one
from dbd4 import parseCSVData
import myPickle
from getPair import mapU2B
def calc_gini(x, use_numpy=True): #follow transformed formula
	"""Return computed Gini coefficient.
	:contact: aisaac AT american.edu
	"""
	xsort = sorted(x) # increasing order
	
	y = np.cumsum(xsort)	
	B = sum(y) / (y[-1] * len(x))
	return 1 + 1./len(x) - 2*B
def getCountDict():
    cnt_dict={}
    for a in AA:
        for b in AA:
            cnt_dict[(a,b)]=0.0
    return cnt_dict
def getVec(d):
    idx=zip(AAcodes,range(len(AAcodes)))
    pp=np.zeros((len(idx),len(d.values()[0])))
    for i,a in enumerate(AAcodes):
        pp[i,:]=d[a]
    return pp  
def getMtx(d):
    idx=zip(AAcodes,range(len(AAcodes)))
    pp=np.zeros((len(idx),len(idx)))
    for i,a in enumerate(AAcodes):
        for j,b in enumerate(AAcodes):
            pp[i,j]=d[(a,b)]+0.1
    return pp      
def setTicksColor(axs):
    for tidx,t in enumerate(axs.get_ticklabels()):
        if tidx<7:
            clr='green'
        elif 7<=tidx<10:
            clr='blue'
        elif 10<=tidx<15:
            clr='orange'
        elif 15<=tidx<18:
            clr='red'
        else:
            clr='black'
        t.set_color(clr)
def addASA(cid,lig,P,A,L0,L1,u2b):
    if lig:            
        lig='l'
    else:            
        lig='r'    
    for r in P:
        skp=0
        try:
            k=three_to_one(L0.R[r].get_resname())
        except KeyError:
            skp=1

            
        try:
            rms=np.max(rmsd[cid+'_'+lig][-2][r])
            
        except KeyError as e:     
            print "Error",e            
            continue
        v0=L0.rASA[r]
        
            
        skp=skp+np.isnan(v0)+np.isnan(u2b[r])+(rms<0)
        
        if skp or np.isnan(L1.rASA[u2b[r]]):
            continue
        v1=L1.rASA[u2b[r]]
        A[k][0]=A[k][0]+1
        A[k][1]=A[k][1]+v0
        A[k][2]=A[k][2]+v1
        A[k][3]=A[k][3]+rms
bdir='../DBD4N/PDBPKL4/'
exfile=bdir+'E_125PN_15_35_50.lbl.pkl'
categs=['All','RB','Med','Hard']#
AAcodes=['Ala','Val','Leu','Ile','Cys','Met','Pro','Phe','Trp','Tyr','Gly','Ser','Thr','Glu','Gln','Arg','Lys','His','Asn','Asp']    
AAcodes=[three_to_one(a.upper()) for a in AAcodes]
#f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN', '1N2C', '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
#f4=['2A5T', '3CPH', '1ZHH', '2ABZ', '1LFD', '2OUL', '1JIW', '2B4J', '1SYX', '1FLE', '1JTG', '2AYO', '4CPA', '1CLV', '1OC0', '1XU1', '1R6Q', '2O3B', '1US7', '3D5S', '1JZD', '1HCF', '1OYV', '2OZA', '1H9D', '2A9K', '2J0T', '2Z0E', '3BP8', '2IDO', '1WDW', '1ZLI', '2VDB', '1RV6', '1FFW', '1F6M', 'BOYV', '1JWH', '2OOR', '1MQ8', '1GL1', '1PVH', '2I9B', '1OFU', '1GXD', '3SGQ', '1JK9', '1ZM4', '1FCC', '2G77', '2J7P', '2FJU']
dbd4=parseCSVData('..\Complete Data\DBD4_data.csv')
fig, axes = plt.subplots(nrows=2, ncols=2)
rmsdfname='rmsd_atomic.mkl' 
rmsd=myPickle.load(rmsdfname)
for subidx,categ in enumerate(categs):

    if categ!='All':
        clist=[cid for cid in dbd4.keys() if dbd4[cid][1].upper()==categ.upper()]
    else:
        clist=dbd4.keys()
    E=getExamplesDBD.loader(exfile)
    Pcnt=getCountDict()
    Ncnt=getCountDict()
    APcnt=dict(zip(AAcodes,[[0,0,0,0] for _ in AAcodes]))
    ANcnt=dict(zip(AAcodes,[[0,0,0,0] for _ in AAcodes]))
    TAC=[]
    
    ofname='propAsabrx_nogly_'+categ+'.prp.mkl'
    if os.path.isfile(ofname) is False:
        for cidx,cid in enumerate(clist):
            print "Currently processing:",categ,cid, round(100*cidx/float(len(clist)),3),'% done.'
            
            Lu=myPDB.loader(bdir+cid+'_l_u.pdb.pkl')
            Lb=myPDB.loader(bdir+cid+'_l_b.pdb.pkl')
            Ru=myPDB.loader(bdir+cid+'_r_u.pdb.pkl')
            Rb=myPDB.loader(bdir+cid+'_r_b.pdb.pkl')
            (Lu2b,_)=mapU2B(Lu.seq,Lu.S2Ri,len(Lu.R),Lb.seq,Lb.S2Ri,len(Lb.R))
            (Ru2b,_)=mapU2B(Ru.seq,Ru.S2Ri,len(Ru.R),Rb.seq,Rb.S2Ri,len(Rb.R))
            
            Pex=E.Pex[cid][0]
            Nex=list(E.getNegEx(cid))
            lPex=np.unique([Pex[i][0] for i in range(len(Pex))])
            rPex=np.unique([Pex[i][1] for i in range(len(Pex))])
            lNex,rNex=E.getNegExSingle(cid)
            #Pex_res=[(three_to_one(Lu.R[Pex[i][0]].get_resname()),three_to_one(Ru.R[Pex[i][1]].get_resname())) for i in range(len(Pex))]
            #Nex_res=[ for i in range(len(Nex))]
            
            rtac=[[Ru.Phi[p],Ru.Psi[p],Rb.Phi[Ru2b[p]],Rb.Psi[Ru2b[p]]] for p in rPex if three_to_one(Ru.R[p].get_resname())!='G']
            ltac=[[Lu.Phi[p],Lu.Psi[p],Lb.Phi[Lu2b[p]],Lb.Psi[Lu2b[p]]] for p in lPex if three_to_one(Lu.R[p].get_resname())!='G']
            TAC.extend(ltac)
            TAC.extend(rtac)
            #pdb.set_trace()
            for a,b in Pex:
                try:
                    k0=(three_to_one(Lu.R[a].get_resname()),three_to_one(Ru.R[b].get_resname()))
                    k1=(three_to_one(Ru.R[b].get_resname()),three_to_one(Lu.R[a].get_resname()))
                except KeyError:    
                    continue           
                Pcnt[k0]=Pcnt[k0]+1
                Pcnt[k1]=Pcnt[k1]+1
            
            for a,b in Nex:
                try:
                    k0=(three_to_one(Lu.R[a].get_resname()),three_to_one(Ru.R[b].get_resname()))
                    k1=(three_to_one(Ru.R[b].get_resname()),three_to_one(Lu.R[a].get_resname()))
                except KeyError:    
                    continue           
                Ncnt[k0]=Ncnt[k0]+1
                Ncnt[k1]=Ncnt[k1]+1
            
            addASA(cid,True,lPex,APcnt,Lu,Lb,Lu2b)
            addASA(cid,False,rPex,APcnt,Ru,Rb,Ru2b)
            addASA(cid,True,lNex,ANcnt,Lu,Lb,Lu2b)
            addASA(cid,False,rNex,ANcnt,Ru,Rb,Ru2b)
            
            #pdb.set_trace()
        myPickle.dump(ofname,(Pcnt,Ncnt,APcnt,ANcnt,TAC))
    else:
        print "Using existing file",ofname
        (Pcnt,Ncnt,APcnt,ANcnt,TAC)=myPickle.load(ofname)
            
    Pm=getMtx(Pcnt)
    Nm=getMtx(Ncnt)
    
    v=np.atleast_2d(np.sum(Nm,axis=0)+np.sum(Pm,axis=0))
    Ex=np.sum(Pm)*((v*v.T)/np.sum((v*v.T)))
    #
    pp=((Pm-Ex)**2)/(Ex)
    
    #pp=(Pm/np.sum(Pm))/(Nm/np.sum(Nm))
    pp=np.log2(Pm/Ex)
    print categ,calc_gini(pp.flatten())
    #pdb.set_trace()
    """
    P=np.sum(Pcnt.values())    
    N=np.sum(Ncnt.values())        
    for k in Pcnt:
        Pcnt[k]=Pcnt[k]/P
        
    for k in Ncnt:
        Ncnt[k]=Ncnt[k]/N    
    
    pdict=dict(zip(Pcnt.keys(),np.array(Pcnt.values())/(1e-10+np.array(Ncnt.values()))))
        
    """
    APm=getVec(APcnt)
    ANm=getVec(ANcnt)
    Ex1=np.sum(APm[:,0])*((ANm[:,0]+APm[:,0])/np.sum(ANm[:,0]+APm[:,0]))
    pm1=APm[:,0]
    pp1=np.log2(pm1/Ex1)
    #pp1=((pm1-Ex1)**2)/(Ex1)
    asax=APm+ANm
    lbl='rASA for all residues'
    if categ!='All':
        asax=APm
        lbl='rASA for binding residues'
    #asax[:,1]=asax[:,1]-asax[:,2]
    asa0=asax[:,1]/asax[:,0]
    asa1=asax[:,2]/asax[:,0]
    if categ=='All':
        vmax=np.max(pp)
        vmin=np.min(pp)
    plt.figure(0)        
    ax=plt.subplot(2,2,subidx+1)
    im=plt.imshow(pp, cmap='jet',interpolation="nearest",vmin=-2.0);#,vmin=vmin,vmax=vmax,,vmin=0.14
    plt.xticks(range(len(AAcodes)),AAcodes);
    plt.yticks(range(len(AAcodes)),AAcodes);
    plt.title(categ)
    plt.colorbar();
    
    setTicksColor(ax.xaxis)
    setTicksColor(ax.yaxis)
    
    plt.figure(1)    
    ax=plt.subplot(2,2,subidx+1)    
    plt.plot(range(len(AAcodes)),pp1,'ro-'); 
    #plt.yscale('log',basex=2)
    
    plt.ylim([-1,+1])
    plt.title(categ)
    plt.ylabel('log(Propesnsity)',color='red')
    
    ax.tick_params(axis='y', colors='red')
    #ax.xaxis.label.set_color('blue')
    plt.grid()
    #if categ=='All':
    ax2 = ax.twinx()
    plt.plot(range(len(AAcodes)),asa0,'bv-',label=lbl+' (unbound)');
    plt.plot(range(len(AAcodes)),asa1,'bs-',label=lbl+' (bound)');
    
    plt.ylim([0,0.6])
    ax2.tick_params(axis='y', colors='blue')
    #ax.xaxis.label.set_color('blue')
    plt.xticks(range(len(AAcodes)),AAcodes);   
    #.ylabel('rASA',color='blue')
    setTicksColor(ax.xaxis)
    plt.legend(loc=4) 
        
    
    plt.figure(2)
    ax=plt.subplot(2,2,subidx+1)
    plt.plot(APm[:,3]/APm[:,0],'ro-',label='Binding residues')
    plt.plot(ANm[:,3]/ANm[:,0],'bs-',label='Non-binding residues')
    plt.ylim([0,3])
    plt.xticks(range(len(AAcodes)),AAcodes); 
    plt.grid()
    plt.title(categ)
    plt.ylabel('RMSD')
    plt.legend(loc=0)
    setTicksColor(ax.xaxis)
#cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
#fig.colorbar(im, cax=cax)
#plt.savefig('../figures/propensities.pdf', format='pdf', dpi=1200)
plt.show()        