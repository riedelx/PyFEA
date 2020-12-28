import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import sys
#sys.path.insert(1, '/libraries')
#import MNcurve as MN
import materials as mat
import utils

def plotFunc(self,x_coordinates,y_coordinates,reverse,title='section layout'):
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    if reverse:y_coordinates=[-i for i in y_coordinates]
    ax.plot(x_coordinates, y_coordinates,'k-')
    for i in self.reinf_sect:
        yTemp=i[2]
        bTemp=self.width(yTemp)
        bSpacing=bTemp/(i[0]+1)
        if reverse:yTemp=-yTemp
        for j in range(i[0]):
            xTemp=-bTemp/2+j*bSpacing+bSpacing
            circ=plt.Circle((xTemp,yTemp), radius=i[1]/2, color='b', fill=False)
            ax.add_patch(circ)
    if reverse:ax.plot(0,-self.centr,'r+',markersize=10,linewidth=8)
    else:ax.plot(0,self.centr,'r+',markersize=10,linewidth=8)
    ax.add_artist(circ)
    ax.set_title(title)
    ax.set_ylabel('Height [mm]')
    ax.set_xlabel('Width [mm]')
    ax.set_aspect('equal', 'box')
    plt.show()

class rcts:
    def __init__(self,Df,Dw,Bf,Bw,reinf_sect):
        self.Dw=Dw
        self.Df=Df
        self.Bw=Bw
        self.Bf=Bf
        reinf_sect = pd.DataFrame(reinf_sect).sort_values(by=[2],ascending=False).astype(int).values.tolist()
        self.reinf_sect=reinf_sect
        self.h=Dw+Df
        area_w = Dw * Bw
        area_f = Df * Bf
        self.area = area_w + area_f
        self.centr = int(self.h-(Df*area_f/2+area_w*(Df+Dw/2))/self.area)
    def plotting(self,title='section layout',reverse=False):
            y_coordinates = [0,self.Dw+0,self.Dw+0,self.h,self.h,self.Dw+0,self.Dw+0,0,0]
            x_coordinates = [self.Bw/2,self.Bw/2,self.Bf/2,self.Bf/2,-self.Bf/2,-self.Bf/2,-self.Bw/2,-self.Bw/2,self.Bw/2]
            plotFunc(self,x_coordinates,y_coordinates,reverse=reverse,title=title)
    def width(self,x):
        if (x >= 0 and x <= self.Dw):
            b=self.Bw
        elif x<=self.h:
            b=self.Bf
        else:
            b=0
        return b

class rcrs:
    def __init__(self, b, d, reinf_sect):
        self.b = b
        self.d = d
        self.h = d
        reinf_sect = pd.DataFrame(reinf_sect).sort_values(by=[2],ascending=False).astype(int).values.tolist()
        self.reinf_sect=reinf_sect
        self.area = b * d
        self.centr = int(d/2)
    def plotting(self,title='section layout',reverse=False):
            y_coordinates = [0,self.d,self.d,0,0]
            x_coordinates = [self.b/2,self.b/2,-self.b/2,-self.b/2,self.b/2]
            plotFunc(self,x_coordinates,y_coordinates,reverse=reverse,title=title)
    def width(self,x):
        if (x >= 0 and x <= self.d):
            b=self.b
        else:
            b=0
        return b
