#============================================
# Common Functions
#============================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import log10, floor

def round_sig(x, sig=2):
    if x ==0: return 0
    else: return round(x, sig-int(floor(log10(abs(x))))-1)

def maxAbs(a,b):
    if np.abs(a) >= np.abs(b):
        return a
    else:
        return b

def matchSign(a,b):
    if b >= 0:
        return a
    else:
        return -a

def find_nearest(array, value): #numpy
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def closestPoint(point, listTo, listFrom):
    index = next(x[0] for x in enumerate(listFrom) if x[1] > point)
    return np.array([float(listTo[index]), float(listFrom[index])])

def df_index(df,val,col_ID): return df.index[df[col_ID] == val].tolist()[0]

def df_value(df,val,col_ID0, col_ID1):
    return df.loc[df_index(df,val,col_ID0)][col_ID1]

def plotBase(grid=True,figsize=(6,4)):
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(111)
    if grid:
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    return fig,ax

# Intersection of two curves
"""
Sukhbinder
5 April 2017
Based on:
"""
def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def intersection(x1,y1,x2,y2):
    """
INTERSECTIONS Intersections of curves.
   Computes the (x,y) locations where two curves intersect.  The curves
   can be broken with NaNs or have vertical segments.
usage:
x,y=intersection(x1,y1,x2,y2)
    Example:
    a, b = 1, 2
    phi = np.linspace(3, 10, 100)
    x1 = a*phi - b*np.sin(phi)
    y1 = a - b*np.cos(phi)
    x2=phi
    y2=np.sin(phi)+2
    x,y=intersection(x1,y1,x2,y2)
    plt.plot(x1,y1,c='r')
    plt.plot(x2,y2,c='g')
    plt.plot(x,y,'*k')
    plt.show()
    """
    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN

    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]

def findExactPoint(curve1, coordinate,limY=True, multiple=False):
    if limY:
        curve2=np.array([[-9E99,np.inf],[coordinate,coordinate]])
    else:
        curve2=np.array([[coordinate,coordinate],[-9E99,np.inf]])
    tupl = (intersection(curve1[0], curve1[1], curve2[0], curve2[1]))
    try:
        return [float(tupl[0]),float(tupl[1])]
    except:
        if multiple:
            return (tupl[0]),(tupl[1])
        else:
            return [float(tupl[0][0]),float(tupl[1][0])]

def centroidX(points,discr=100):
    #points=[[x1,y1],[x2,y2]]
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    points=np.array([x,y])
    area=np.trapz([y[0],y[-1]], x=[x[0],x[-1]])
    L=x[-1]-x[0]
    L_incr=L/discr
    moment=0
    for i in range(discr):
        if i == discr-1:
            x2=x[-1]
        else:
            x2=(i+1)*L_incr
        y2=findExactPoint(points, x2,limY=False)[1]
        x1=(i)*L_incr
        y1=findExactPoint(points, x1,limY=False)[1]
        area_temp=np.trapz([y1,y2], x=[x1,x2])
        moment+=area_temp*(x2-L_incr/2)
    return moment/area
