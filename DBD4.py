# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 23:55:30 2013
Helper functions for Docking benchmark data set
Returns information about the DBD4.py originally in its excel file ..\Complete Data\DBD4_data.csv
@author: root
"""
import csv
import urllib2
import myPickle
def parseCSVData(dbd4file):
    """
    Returns a dictionary object for the file ..\Complete Data\DBD4_data.csv
        cid  : complete cid, category, lid, rid, rmsd, buried_asa
    """
    dbd4={}
    nhdr=1
    i=0    
    with open(dbd4file, 'rb') as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            if len(row)>0:
                if len(row[0]):                
                    if row[0][0]=='#':
                        categ=row[0][1:]
                        continue
                    elif len(row)==5:
                        i=i+1
                        if i>nhdr:
                            cid,lid,rid,rmsd,asa=row                        
                            dbd4[cid[:4]]=(cid,categ,lid,rid,float(rmsd),float(asa))    
    return dbd4                     
    


def getStoichiometry(cid):   
    
    website = urllib2.urlopen('http://www.rcsb.org/pdb/explore/explore.do?structureId='+cid)
    wtxt = website.read()
    #website_text = retrieveWebPage('http://www.rcsb.org/pdb/explore/explore.do?structureId=+'+cid).read()
    
    idx=wtxt.find('Stoichiometry:')    
    stxt=wtxt[idx:idx+100]
    stch=stxt[(stxt.find('<b>')+3):(stxt.find('</b>'))]
    
    idx=wtxt.find('Symmetry:')    
    stxt=wtxt[idx:idx+100]
    symm=stxt[(stxt.find('<b>')+3):(stxt.find('</b>'))]
    return stch,symm

def kMer(stch):
    if len(stch):
        return int(stch[stch.find('-mer')-2:][0:2])
    else:
        return 0
        
if __name__=="__main__":
    
    
    dbd4file='..\..\Complete Data\DBD4_data.csv'
    dbd4=parseCSVData(dbd4file)
    dbd4s=myPickle.load('../../Complete Data/DBD4_data_stocihiometry.mkl')
    from getExamplesDBD_breakup import *
    E=getExamplesDBD.loader('../../DBD4CSPKL/PKL/ENS_15_35_50.lbl.pkl')
    dbd4es={}
    for cid in E.Pex:
        scid=cid[0][:4]
        dbd4es[cid]=dbd4s[scid]+(E.Pex[cid][1],E.Pex[cid][2],E.Pex[cid][1]*E.Pex[cid][2],len(E.Pex[cid][0]),kMer(dbd4s[scid][6]))
    with open('DBD4_broken.csv', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for cid in dbd4es:
            spamwriter.writerow([cid]+list(dbd4es[cid]))  
#    bdir='../../Complete Data/CPDB/'
#    from BISEPutils import fetchPDB
#    for c in dbd4:
#        try:
#            fetchPDB(c,bdir=bdir)
#            print "Done",c
#        except Exception as e:
#            print "Error fetching",c,e
#    for cid in dbd4.keys():
#        try:
#            stx=getStoichiometry(cid)
#            print cid,stx
#            dbd4[cid]+=stx
#        except:
#            print "Failed",cid
#    myPickle.dump('DBD4_data_stocihiometry.mkl',dbd4)
#    #with open('..\Complete Data\DBD4_data_stocihiometry.csv','wb') as f:
#    #    w = csv.writer(f)
#    #    w.writerows(dbd4.items())

#    dbd4=myPickle.load('../../Complete Data/DBD4_data_stocihiometry.mkl')
#    R=myPickle.load('SGD_DBD4.res.pkl')
#    J={}
#    for cid in dbd4:
#        try:
#            stch=dbd4[cid][-2]            
#            J[cid]=dbd4[cid]+(kMer(stch),)+R[cid][0][2:]+R[cid][0][1]+R[cid][0][0]
#        except:
#            continue
#    
#
#    with open('SGD_DBD4.csv', 'wb') as csvfile:
#        spamwriter = csv.writer(csvfile, delimiter=',')
#        for cid in J:
#            spamwriter.writerow([cid]+list(J[cid]))    
        
        
    # CONSTRUCT HISTOGRAM OF k-mers
    
    