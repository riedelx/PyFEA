import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import sys
#sys.path.insert(1, '/libraries')
#import MNcurve as MN
import materials as mat
import utils

def plotFunc(self,x_coordinates,y_coordinates):
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    ax.plot(x_coordinates, y_coordinates,'k-')
    #print(x_coordinates)
    #print(y_coordinates)
    for i in self.reinf_sect:
        #print(i)
        xTemp=i[2]
        bTemp=self.width(xTemp)
        bSpacing=bTemp/(i[0]+1)
        for j in range(i[0]):
            yTemp=-bTemp/2+j*bSpacing+bSpacing
            #print(j,bTemp,xTemp,yTemp,i[1])
            circ=plt.Circle((xTemp,yTemp), radius=i[1]/2, color='b', fill=False)
            ax.add_patch(circ)
    ax.plot(self.centr,0,'r+',markersize=10,linewidth=8)
    ax.add_artist(circ)
    ax.set_title('section layout')
    ax.set_aspect('equal', 'box')
    plt.show()

class rcts:
    def __init__(self,Df,Dw,Bf,Bw,reinf_sect,plotting=True):
        self.Dw=Dw
        self.Df=Df
        self.Bw=Bw
        self.Bf=Bf
        self.reinf_sect=reinf_sect
        self.h=Dw+Df
        area_w = Dw * Bw
        area_f = Df * Bf
        self.area = area_w + area_f
        self.centr = int(self.h-(Df*area_f/2+area_w*(Df+Dw/2))/self.area)
        if plotting:
            x_coordinates = [0,self.Dw+0,self.Dw+0,self.h,self.h,self.Dw+0,self.Dw+0,0,0]
            y_coordinates = [self.Bw/2,self.Bw/2,self.Bf/2,self.Bf/2,-self.Bf/2,-self.Bf/2,-self.Bw/2,-self.Bw/2,self.Bw/2]
            plotFunc(self,x_coordinates,y_coordinates)
    def width(self,x):
        if (x >= 0 and x <= self.Dw):
            b=self.Bw
        elif x<=self.h:
            b=self.Bf
        else:
            b=0
        return b

class rss:
    def __init__(self, b, d, reinf_sect, plotting = True):
        self.b = b
        self.d = d
        self.h = d
        self.reinf_sect=reinf_sect
        self.area = b * d
        self.centr = int(d/2)
        if plotting:
            x_coordinates = [0,d,d,0,0]
            y_coordinates = [b/2,b/2,-b/2,-b/2,b/2]
            plotFunc(self,x_coordinates,y_coordinates)
    def width(self,x):
        if (x >= 0 and x <= self.d):
            b=self.b
        else:
            b=0
        return b
