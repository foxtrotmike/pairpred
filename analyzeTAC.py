# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 22:44:32 2013
Plots the binding associated changes in torsion angles after clustering analysis for a number of proteins
@author: root
"""
import numpy as np
import myPickle
from scipy.cluster.vq import *
import matplotlib
import matplotlib.pyplot as plt
#from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay

def voronoi(P):
    delauny = Delaunay(P)
    triangles = delauny.points[delauny.vertices]

    lines = []

    # Triangle vertices
    A = triangles[:, 0]
    B = triangles[:, 1]
    C = triangles[:, 2]
    lines.extend(zip(A, B))
    lines.extend(zip(B, C))
    lines.extend(zip(C, A))
    lines = matplotlib.collections.LineCollection(lines, color='r')
    #plt.gca().add_collection(lines)

    circum_centers = np.array([triangle_csc(tri) for tri in triangles])

    segments = []
    for i, triangle in enumerate(triangles):
        circum_center = circum_centers[i]
        for j, neighbor in enumerate(delauny.neighbors[i]):
            if neighbor != -1:
                segments.append((circum_center, circum_centers[neighbor]))
            else:
                ps = triangle[(j+1)%3] - triangle[(j-1)%3]
                ps = np.array((ps[1], -ps[0]))

                middle = (triangle[(j+1)%3] + triangle[(j-1)%3]) * 0.5
                di = middle - triangle[j]

                ps /= np.linalg.norm(ps)
                di /= np.linalg.norm(di)

                if np.dot(di, ps) < 0.0:
                    ps *= -1000.0
                else:
                    ps *= 1000.0
                segments.append((circum_center, circum_center + ps))
    return segments

def triangle_csc(pts):
    rows, cols = pts.shape

    A = np.bmat([[2 * np.dot(pts, pts.T), np.ones((rows, 1))],
                 [np.ones((1, rows)), np.zeros((1, 1))]])

    b = np.hstack((np.sum(pts * pts, axis=1), np.ones((1))))
    x = np.linalg.solve(A,b)
    bary_coords = x[:-1]
    return np.sum(pts * np.tile(bary_coords.reshape((pts.shape[0], 1)), (1, pts.shape[1])), axis=0)
    
    
categs=['All','Hard']
for pidx,categ in enumerate(categs):
    fname='Data_out/propAsabrx_nogly_'+categ+'.prp.mkl'
    (Pcnt,Ncnt,APcnt,ANcnt,TAC)=myPickle.load(fname)
    TAC=np.array(TAC)
    TAC=TAC[~np.any(TAC>180,axis=1),:]
    if categ=='All':
        Nc=60
        niter=2000
        res0, _ = kmeans2(np.vstack((TAC[:,:2],TAC[:,2:])),Nc,iter=niter,minit='points')
        
        res=np.zeros((Nc**2,4))
        k=0
        for i in range(Nc):
            for j in range(Nc):
                res[k,:]=np.hstack((res0[i,:],res0[j,:]))
                k=k+1
    idx = vq(TAC, res)[0]

    cnt=dict(zip(range(res.shape[0]),[0 for _ in range(res.shape[0])]))
    for i in idx:
        cnt[i]=cnt[i]+1.0
    N=np.sum(cnt.values())        
    plt.figure(pidx)#plt.subplot(2,2,pidx+1)   
    ax=plt.gca()
    plt.plot(TAC[:,0],TAC[:,1],'b.',markersize=1.0)  
    plt.plot(TAC[:,2],TAC[:,3],'r.',markersize=1.0)  
    for i in range(res.shape[0]):
        dphi=res[i,2]-res[i,0]
        dpsi=res[i,3]-res[i,1]   
        plt.plot(res[i,0],res[i,1],'bo')
        if ((np.abs((dphi + 180) % 360 - 180)+np.abs((dpsi + 180) % 360 - 180)) > 60) and cnt[i]>=2 :#and not np.any((res[i,:]>180) + (res[i,:]<-180))
            
            #plt.plot(res[i,2],res[i,3],'ro')
            #lw=np.min((8,np.exp(500*cnt[i]/N)))
            if cnt[i]>=0 and cnt[i]<4:
                clr='0.5'
                lw=0.5
            elif cnt[i]>=4 and cnt[i]<8:
                clr='g'
                lw=2.0
            else:
                clr='k'
                lw=4.0
            #print res[i,:],cnt[i]
           # plt.arrow( res[i,0], res[i,1], dphi,dpsi ,fc="k", ec="g",head_width=5, head_length=10,linewidth=200*cnt[i]/N)
            ax.annotate("",
            xy=(res[i,0], res[i,1]), xycoords='data',
            xytext=(res[i,2],res[i,3]), textcoords='data',
            arrowprops=dict(arrowstyle="->", #linestyle="dashed",
                            color=clr,linewidth=lw,
                            patchB=None,
                            shrinkB=0,
                            connectionstyle="arc3,rad=0.3",
                            ),
            )
            #matplotlib.patches.FancyArrowPatch.set_connectionstyle("arc,angleA=0,armA=30,rad=10")
            #plt.plot(res[i,0],res[i,1],'bo')
            #plt.plot(res[i,2],res[i,3],'ro')
    #pdb.set_trace()
    
    plt.xlim([-180,180])    
    plt.ylim([-180,180])
    #plt.axis('equal')
    plt.grid()
    plt.title(categ)
    plt.xlabel('$\Phi$')
    plt.ylabel('$\Psi$')
    segments = voronoi(res0)
    lines = matplotlib.collections.LineCollection(segments, color='k')
    ax.add_collection(lines)
    #segments = voronoi(res[:,2:])
    #lines = matplotlib.collections.LineCollection(segments, color='r')
    ax.add_collection(lines)
    
    #voronoi_plot_2d(res)
plt.show()    