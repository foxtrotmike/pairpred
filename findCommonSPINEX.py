# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:29:53 2013
Find the proteins that are common between SPINEX training set and DBD 4.0
@author: root
"""
import csv
# import pdb
lstfile='..\Tools\spineXpublic\list.2640'
dbd4file='..\Complete Data\DBD4_data.csv'
lst=[]
with open(lstfile, 'r') as f:
    for ln in f:
        lns=ln.strip()
        pid=lns[:4]+'_'+lns[4:]
        lst.append(pid)

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

def catX(f):
    if len(f[4:])>1:
        pid=f[:4]
        rlist=[]
        for c in f[5:]:
            rlist.append(pid+'_'+c)
        return rlist            
    else:
        return [f]                        
f3=['1SBB', '1JPS', '2HMI', '1GHQ', '1KTZ', '1K74', '1D6R', '2SIC', '1GPW', '1XD3', '1EAW', '1VFB', '7CEI', '1E4K', '1I4D', '1H1V', '2PCC', '1FQ1', '2HLE', '1FQJ', '1S1Q', '2OOB', '1UDI', '1KLU', '1WQ1', '1CGI', '1ATN',  '1GP2', '1FAK', '1NW9', '1GLA', '1GRN', '2HRK', '1AZS', '1JMO', '1PXV', '1EWY', '1RLB', '1DQJ', '2BTF', '2I25', '1I2M', '1BUH', '1BGX', '1ML0', '1EFN', '1DFJ', '1Y64', '2UUY', '1MAH', '1BVK', '1BVN', '1EER', '1MLC', '1NSN', '1AK4', '1A2K', '1QFW', '2H7V', '1T6B', '1KAC', '1YVB', '1J2J', '1QA9', '1AHW', '2OT3', '2FD6', '2AJF', '1K4C', '1NCA', '1OPH', '1XQS', '1B6C', '1PPE', '2O8V', '1HIA', '1Z0K', '1R0R', '1WEJ', '1ACB', '1KXP', '1KXQ', '1R8S', '1IRA', '1GCQ', '1F51', '2B42', '2HQS', '1AKJ', '2JEL', '1KKL', '1FC2', '1E96', '1N8O', '2MTA', '2VIS', '1IB1', '1E6J', '1Z5Y', '1EZU', '1TMQ', '2C0L', '1E6E', '1IQD', '1ZHI', '1M10', '2NZ8', '1AY7', '1HE8', '1IJK', '1HE1', '1FSK', '1F34', '2SNI', '1BJ1', '2CFH', '1BKD', '1DE4', '1IBR', '1I9R', '1K5D', '1AVX']
prot4=[]
for cid in dbd4.keys():
    a=[]
    a.extend(catX(dbd4[cid][2]))
    a.extend(catX(dbd4[cid][3]))
    ts=0
    for b in a:
        if b in lst:
            ts=1
    if ts:
        continue            
    prot4.append(cid)

print list(set(f3).intersection(prot4))    