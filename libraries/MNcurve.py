import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import sys
#sys.path.insert(1, '/libraries')
import sections as sect
import materials as mat
import utils


def findPointZero(curve,epsilon,limY=False):
    try: return utils.findExactPoint(curve, epsilon,limY=limY)[1]
    except: return 0

class rctsSect:
    def __init__(self,Df,Dw,Bf,Bw,reinf_sect):
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
        #self.xmin = -int(self.h-(Df*area_f/2+area_w*(Df+Dw/2))/self.area)
        #self.xmax = int((Df*area_f/2+area_w*(Df+Dw/2))/self.area)
    def width(self,x):
        if (x >= 0 and x <= 500):
            b=self.Bw
        elif x<=self.h:
            b=self.Bf
        else:
            b=0
        return b
    def plot(self):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x_coordinates = [0,self.Dw+0,self.Dw+0,self.h,self.h,self.Dw+0,self.Dw+0,0,0]
        y_coordinates = [self.Bw/2,self.Bw/2,self.Bf/2,self.Bf/2,-self.Bf/2,-self.Bf/2,-self.Bw/2,-self.Bw/2,self.Bw/2]
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

class MNclass:
    def __init__(self,concr,reinf,section):
        self.h=section.h
        self.section=section
        self.concr=concr
        self.reinf=reinf
        stl1_df = reinf.stress_df()
        self.stl1_df=pd.concat([-stl1_df.iloc[::-1],stl1_df])[::-1]

    def materials(self,plotting=True):
        self.stl1=self.plotStress(self.stl1_df,lbl="stl1",title="Reinforcement",xlim=(None,None),ylim=(None,None),plotting=plotting)
        self.stl1=np.array(self.stl1)
        self.con1_df=pd.DataFrame([[0,self.concr.epsilon_2t],[self.concr.ft,self.concr.epsilon_1t],
                              [0,0],[-self.concr.fc1,-self.concr.epsilon_1c],[-self.concr.fc2,-self.concr.epsilon_2c]],columns=['stress','strain'])
        self.con1=self.plotStress(self.con1_df,lbl="con1",title="C50/60")
        self.con1=np.array(self.con1)

    def plotStress(self,curve,title="",lbl="",xlim=(None,None),ylim=(None,None),plotting=True):
        if plotting:
            fig,ax = utils.plotBase()
            ax.plot(curve['strain'],curve['stress'],'-', linewidth=2, markersize=5,label=lbl)
            ax.legend(loc='lower right')
            ax.set_title(title)
            ax.set_xlabel('Strain')
            ax.set_ylabel('Stress [MPa]')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            plt.show()
        return curve['strain'],curve['stress']
    def epsilonBuildX(self,x_NA,hogging,eps_u=0.0035,plotting=False):
        if hogging:
            epsH=eps_u
            eps0=eps_u*(1-h/(h-x_NA))
        else:
            eps0=eps_u
            epsH=eps_u*(1-h/(x_NA))
        if plotting:
                fig,ax = utils.plotBase()
                ax.plot([0,h],[eps0,epsH],'-', linewidth=2, markersize=5)
                ax.plot([0,h],[0,0],'-', linewidth=2, markersize=5,color='black')
                ax.set_title('Strain distribution within section')
                ax.set_xlabel('Distance [mm]')
                ax.set_ylabel('Strain []')
                ax.set_xlim(0,h)
                ax.set_ylim(None,None)
                plt.show()
        return np.array([[0,h],[eps0,epsH]])
    def epsilonBuildEps(self,epsH,eps0,plotting=False):
        if plotting:
                fig,ax = utils.plotBase()
                ax.plot([0,self.h],[eps0,epsH],'-', linewidth=2, markersize=5)
                ax.plot([0,self.h],[0,0],'-', linewidth=2, markersize=5,color='black')
                ax.set_title('Strain distribution within section')
                ax.set_xlabel('Distance [mm]')
                ax.set_ylabel('Strain []')
                ax.set_xlim(0,self.h)
                ax.set_ylim(None,None)
                plt.show()
        return np.array([[0,self.h],[eps0,epsH]])
    def epsilonFunc(self,x,epsilon): return utils.findExactPoint(epsilon, x,limY=False)[1]
    def calc(self,eps0,epsH,plotting=False,n_layers=100,eps_u=0.0035):
        epsilon=self.epsilonBuildEps(eps0=eps0,epsH=epsH,plotting=plotting)
        #strain_conLim=.0035
        h_i=self.h/n_layers

        # Steel
        f_s=[]
        m_s=0
        sigma_s=[]
        x_s=np.array(self.section.reinf_sect).T[2]
        for i in self.section.reinf_sect:
            eps_s=self.epsilonFunc(i[2],epsilon)
            sigma= findPointZero(self.stl1,eps_s)
            sigma_s.append(sigma)
            A=np.pi*i[1]**2/4*i[0]
            f_s_i=sigma*A
            f_s.append(f_s_i)
            m_s+=f_s_i*i[2]
        # Steel stress distribution
        sigma_sEnv=[]
        x_sEnv=[]
        for i in range(n_layers):
            x_i=i*h_i
            x_sEnv.append(x_i)
            b_i=self.section.width(x_i)
            e=self.epsilonFunc(x_i,epsilon)
            s=findPointZero(self.stl1,e)
            sigma_sEnv.append(s)
        if plotting:
            fig,ax = utils.plotBase()
            ax.bar(x_s,sigma_s,width =5)
            ax.plot(x_sEnv,sigma_sEnv,'--', linewidth=2, markersize=5)
            ax.plot([0,self.h],[0,0],'-', linewidth=2, markersize=5,color='black')
            ax.set_title('Steel stress distribution within section')
            ax.set_xlabel('Distance [mm]')
            ax.set_ylabel('Stress [MPa]')
            ax.set_xlim(0,self.h)
            ax.set_ylim(None,None)
            plt.show()

        # Concrete
        sigma_con=[]
        x_con=[]
        f_con=[]
        m_con=0
        for i in range(n_layers):
            x_i=i*h_i
            x_con.append(x_i)
            b_i=self.section.width(x_i)
            e=self.epsilonFunc(x_i,epsilon)
            s=findPointZero(self.con1,e)
            sigma_con.append(s)
            f_con_i=s*b_i*h_i
            f_con.append(f_con_i)
            m_con+=f_con_i*x_i
            #print('e: {0}, s: {1}'.format(e,s))
        if plotting:
            fig,ax = utils.plotBase()
            ax.plot(x_con,sigma_con,'-', linewidth=2, markersize=5)
            ax.plot([0,self.h],[0,0],'-', linewidth=2, markersize=5,color='black')
            ax.set_title('Concrete stress distribution within section')
            ax.set_xlabel('Distance [mm]')
            ax.set_ylabel('Stress [MPa]')
            ax.set_xlim(0,self.h)
            ax.set_ylim(None,None)
            plt.show()
        f_tot=sum(f_con)+sum(f_s)
        m_tot=m_con+m_s
        return f_tot,m_tot,f_s,f_con#,f_s,m_con,m_s#sum(f_con)+sum(f_s),m_con+m_s,sum(f_con),sum(f_s),m_con,m_s
