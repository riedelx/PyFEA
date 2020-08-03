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

class MNclass:
    def __init__(self,concr,reinf,section):
        self.h=section.h
        self.section=section
        self.concr=concr
        self.reinf=reinf

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
    def calcX0(self,eps0,x_NA,plotting=False,n_layers=100):
        if x_NA == 0:
            epsH = 0
        elif x_NA>self.h:
            eps0 = eps0/2*x_NA/(x_NA-0.5*self.h)
            epsH = - eps0 * (self.h - x_NA) / x_NA
        else:
            epsH = - eps0 * (self.h - x_NA) / x_NA
        return self.calc(eps0,epsH,plotting=plotting,n_layers=n_layers)

    def calcXH(self,epsH,x_NA,plotting=False,n_layers=100):
        if x_NA == 0:
            eps0 = 0
        elif x_NA>self.h:
            epsH = epsH/2*x_NA/(x_NA-0.5*self.h)
            eps0 = - epsH * (self.h - x_NA) / x_NA
        else:
            eps0 = - epsH * (self.h - x_NA) / x_NA
        return self.calc(eps0,epsH,plotting=plotting,n_layers=n_layers)
    def calc(self,eps0,epsH,plotting=False,n_layers=100):
        epsilon=self.epsilonBuildEps(eps0=eps0,epsH=epsH,plotting=plotting)
        #strain_conLim=.0035
        h_i=self.h/n_layers

        # Steel
        f_s=[]
        m_s=0
        sigma_s=[]
        eps_s=[]
        x_s=np.array(self.section.reinf_sect).T[2]
        for i in self.section.reinf_sect:
            eps=self.epsilonFunc(i[2],epsilon)
            eps_s.append(eps)
            sigma= findPointZero(self.reinf.np,eps)
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
            s=findPointZero(self.reinf.np,e)
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
            s=findPointZero(self.concr.np,e)
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
        m_tot=m_con+m_s - 0.5*f_tot*self.h
        if plotting:
            print('Total axial force: {} kN'.format(int(f_tot/1E3)))
            print('Total moment: {} kNm'.format(int(m_tot/1E6)))
        return f_tot,m_tot,f_s,f_con,eps_s,sigma_s
    def mnCurve(self,xRatio=[0.16,0.2,0.3,0.4,0.5,0.8,0.9,1,1E99],n_layers=100,epsU=-0.0035,reverseMoment=False):
        F=[]
        M=[]
        xRatio=[i*self.h for i in xRatio]
        f_tot,m_tot,f_s,f_con,eps_s,sigma_s=self.calcX0(eps0=self.reinf.epsilon_u,x_NA=1E99, plotting=False, n_layers=n_layers)
        F.append(int(f_tot/1E3))
        M.append(int(m_tot/1E6))
        for i in xRatio:
            #print(i)
            f_tot,m_tot,f_s,f_con,eps_s,sigma_s=self.calcX0(eps0=epsU,x_NA=i, plotting=False, n_layers=n_layers)
            F.append(int(f_tot/1E3))
            M.append(int(m_tot/1E6))
        #for i in xRatio[:-1]:
        for i in xRatio[-2::-1]:
            #print(-i)
            f_tot,m_tot,f_s,f_con,eps_s,sigma_s=self.calcXH(epsH=epsU,x_NA=i, plotting=False, n_layers=n_layers)
            F.append(int(f_tot/1E3))
            M.append(int(m_tot/1E6))
        f_tot,m_tot,f_s,f_con,eps_s,sigma_s=self.calcX0(eps0=self.reinf.epsilon_u,x_NA=1E99, plotting=False, n_layers=n_layers)
        F.append(int(f_tot/1E3))
        M.append(int(m_tot/1E6))
        mnInteraction = pd.DataFrame(np.array([F,M]).T,columns=['F','M'])#.sort_values(by=['x'])
        if reverseMoment:
            mnInteraction['M']=-mnInteraction['M']
        fig,ax = utils.plotBase()
        ax.plot(mnInteraction['M'],mnInteraction['F'],'-o', linewidth=2, markersize=5)
        ax.set_title('M-N interaction diagram')
        ax.set_xlabel('Moment [kNm]')
        ax.set_ylabel('Axial load [kN]')
        plt.show()
        return mnInteraction
