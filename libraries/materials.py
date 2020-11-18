import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import utils

def plotStress(self,curve,title="",lbl="",xlim=(None,None),ylim=(None,None),plotting=True,legend=True):
    if plotting:
        fig,ax = utils.plotBase()
        ax.plot(curve['strain'],curve['stress'],'-', linewidth=2, markersize=5,label=lbl)
        if legend:
            ax.legend(loc='lower right')
        ax.set_title(title)
        ax.set_xlabel('Strain')
        ax.set_ylabel('Stress [MPa]')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()
    return curve['strain'],curve['stress']

class con1:
    def __init__(self, ID, fc1, length, epsilon_t2 = 0.001, fc2_factor = 0.1, ft_factor = 1, characteristic = True,plotting=True,title="con1",Ec2 = '',Et2 = '',strain_prec=5,legend=True):
        self.resid_str = fc2_factor
        self.ID = ID
        self.fc1 = round(fc1, 1)
        self.fc2 = round(self.resid_str * fc1, 1)
        self.length = length
        if characteristic:
            self.fcm = round(fc1+8,2)
        else:
            self.fcm = fc1
        if fc1 <= 50:
            self.ft = round(ft_factor * 0.3 * self.fcm ** (2/3), 1)
        else:
            self.ft = round(ft_factor * 2.12*np.log(1+0.1*self.fcm), 1)
        self.Gf = round(73 * self.fcm**0.18/1000, 3)
        self.Ec0 = round(int(21500*(self.fcm/10)**(1/3)),-2)
        self.poisson = 0.2
        self.Gc = round(250 * self.Gf, 1)
        self.epsilon_1c = round(5 * fc1 / self.Ec0 /3, strain_prec)
        self.Ec1 = int(round(fc1 / self.epsilon_1c, -2))
        self.alpha = min(max(0,round((self.Ec0 - self.Ec1)/self.Ec1,2)),1)
        if Ec2 != '':
            self.Ec2 = Ec2
            self.epsilon_2c = round((self.fc1-self.fc2)/-Ec2+self.epsilon_1c, strain_prec)
        else:
            self.epsilon_2c = round(self.epsilon_1c + 3 * self.Gc / (2 * length * fc1), strain_prec)
            self.Ec2 = - int(round((1 - self.resid_str) * fc1 /(self.epsilon_2c - self.epsilon_1c), -2))
        self.Et1 = self.Ec0
        self.epsilon_1t = round(self.ft / self.Et1, strain_prec)
        self.epsilon_2t = epsilon_t2
        if Et2 != '':
            self.Et2 = Et2
            self.epsilon_2t = round((self.ft)/-Et2+self.epsilon_1t, strain_prec)
        else:
            self.epsilon_2t = epsilon_t2
            self.Et2 = - int(round(self.ft /(self.epsilon_2t - self.epsilon_1t), -2))

        self.df=pd.DataFrame([[0,self.epsilon_2t],[self.ft,self.epsilon_1t],
                              [0,0],[-self.fc1,-self.epsilon_1c],[-self.fc2,-self.epsilon_2c]],columns=['stress','strain'])
        self.np=plotStress(self,self.df,lbl="con1",plotting=plotting,title=title,legend=legend)
        self.np=np.array(self.np)

    def adaptic_print(self):
        line = utils.str_joint([self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,
                                      self.ft, self.Et2, self.alpha])
        return line

    def data_frame(self):
        data = np.array([[self.ID, self.length, self.fc1, self.fc2, self.ft, self.Ec0, self.Ec1,
                          self.Ec2, self.Et1, self.Et2, self.Gf, self.Gc, self.epsilon_1c,
                          self.epsilon_2c, self.epsilon_1t, self.epsilon_2t,  self.alpha]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
                      '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
                      '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
                      '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha$$']
        return df

class stl1:
    def __init__(self, ID, E1, fy, fu, epsilon_u,plotting=True,title="stl1",tension=True,compression=True,legend=True):
        self.ID = ID
        self.E1 = E1
        self.fy = fy
        self.fu = fu
        self.epsilon_u = (epsilon_u)
        self.epsilon_y = round(fy / E1,5)
        self.E2 = round((fu - fy) / (epsilon_u - self.epsilon_y),1)
        self.mu = round(self.E2 / E1,7)

        self.df = pd.DataFrame([[0,0],[self.epsilon_y,self.fy],[self.epsilon_u,self.fu]],columns=['strain','stress'])
        if tension and compression:
            self.df=pd.concat([-self.df.iloc[::-1],self.df])[::-1]
        elif compression:
            self.df=-self.df

        self.np=plotStress(self,self.df,lbl="stl1",title=title,xlim=(None,None),ylim=(None,None),plotting=plotting,legend=legend)
        self.np=np.array(self.np)

    def data_frame(self):
        data = np.array([[self.ID, self.E1, self.E2, self.fy, self.fu, self.epsilon_y, self.epsilon_u,
                          self.mu]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$E_{1}[MPa]$$', '$$E_{2}[MPa]$$', '$$f_{y}[MPa]$$', '$$f_{u}[MPa]$$',
                      '$$e_{y}$$','$$e_{u}$$','$$mu$$']
        return df

class esb1:
    def __init__(self, ID, fu, epsilon_u=0.0035, plotting=True,title="esb1"):
        self.ID = ID
        self.fu = fu
        self.epsilon_u2 = round(epsilon_u,4)
        self.epsilon_u1 = 0.2*self.epsilon_u2

        self.df = pd.DataFrame([[0,0],[self.epsilon_u1,0],[self.epsilon_u1,self.fu],[self.epsilon_u2,self.fu]],columns=['strain','stress'])
        self.df=-self.df

        self.np=plotStress(self,self.df,lbl="stl1",title=title,xlim=(None,None),ylim=(None,None),plotting=plotting)
        self.np=np.array(self.np)

    def data_frame(self):
        data = np.array([[self.ID, self.fu, self.epsilon_u1, self.epsilon_u2]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$f_{u1}[MPa]$$', '$$e_{u1}$$', '$$e_{u2}$$']

class stmdl2:
    # this is con1 ADAPTIC model
    def __init__(self,ec0,muec1,strnc1,stresc1,et0,muet1,strnt1,pset,crks,quad):#pseto,crkso,
        # This subroutine calculates the stress at a monitoring point for
        # material MODEL(2).

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain
        self.ec0=ec0 # Secant compressive stiffness
        self.muec1=muec1 # Compressive softening stiffness
        self.strnc1=strnc1 # strain at residual compressive strength
        self.stresc1=stresc1 # residual compressive strength
        self.et0=et0 # Tensile stiffness
        self.muet1=muet1 # Tensile softening stiffness
        self.strnt1=strnt1 # strain at peak tensile strength ?
        self.pseto=pseto # plastic strain in compression at the start of the step,
                         # represents the intersection of the unloading branch with the strain axis
        self.crkso=crkso # plastic strain in tension at the start of the step
        self.pset=pseto
        self.crks=crkso
        self.quad=quad

    def stress(self, strn):
        self.strn=strn
        if(self.strn<=self.pseto):      # Compressive strain increment
            # NOTE: quad is the relative difference between the initial
            #       compressive tangent modulus and the secant modulus

            self.ec0t=(1+self.quad)*self.ec0 # initial tangent modulus in compression

            if(self.quad>0): # implies quadratic initial compressive response
                self.strnc0=(self.stresc1-self.muec1*self.strnc1)/(self.ec0-self.muec1)

            # Obtain the stress assuming elastic conditions, and determine
            # the force from the limiting curve
            self.strese=self.ec0t*(self.strn-self.pseto) # elastic stress based on ec0t
            # this gives stress in intial compressive response if quad = 0

            if(self.quad>0 and self.strn>self.strnc0): # initial quadratic response
                # pseto>=strn>strnc0

                #quadratic formulation for stress:
                self.stresl=self.strn*(self.ec0t+(self.ec0-self.ec0t)*self.strn/self.strnc0)
                self.etan=self.ec0t+2*(self.ec0-self.ec0t)*self.strn/self.strnc0

            elif(strn>strnc1): # softening branch

                self.stresl=self.stresc1+self.muec1*(self.strn-self.strnc1)
                self.etan=self.muec1

            else: # residual compressive strength

                self.stresl=self.stresc1
                self.etan=0.0

            # Establish the stress and the plastic set

            if(self.strese>self.stresl):

                self.stres=self.strese
                self.pset=self.pseto
                self.etan=self.ec0t

            else:

                self.stres=self.stresl
                self.pset=self.strn-self.stresl/self.ec0t

            self.crks=self.crkso

        elif(self.et0==0.0 or self.strn<self.crkso+self.pseto):   # Cracked zone
            self.stres=0.0
            self.pset=self.pseto
            self.crks=self.crkso
            self.etan=0

        else:   # Tensile strain increment

            # Obtain the stress assuming elastic conditions, and
            # determine the force from the limiting curve

            self.strese=self.et0*(self.strn-(self.pseto+self.crkso))

            self.stresl=self.muet1*(self.strn-(self.pseto+self.strnt1))+self.et0*self.strnt1 # my modification

            if(self.stresl>0.0):

                self.etan=self.muet1

            else:

                self.stresl=0.0
                self.etan=0.0

            # Establish the stress and the plastic set

            if(self.strese<self.stresl): # initial elastic tensile response

                self.stres=self.strese
                self.etan=self.et0
                self.crks=self.crkso

            else: # tensile softening

                self.stres=self.stresl
                self.crks=self.strn-self.pseto-self.stresl/self.et0

            self.pset=self.pseto
        return self.stres

    def plot(self,title='stmdl2',step_size=0.001,legend=True,ranges=''):
        if ranges=='':ranges=[-np.absolute(self.strnc1),np.absolute(self.strnt1)]
        strain=np.arange(ranges[0],ranges[1],step_size)
        stress=np.array([self.stress(i) for i in strain])
        strain=strain.reshape(len(strain),1)
        stress=stress.reshape(len(stress),1)
        self.df=pd.DataFrame(np.hstack((stress,strain)),columns=['stress','strain'])
        self.np=plotStress(self,self.df,lbl="stmdl2",plotting='True',title=title,legend=legend)
        self.np=np.array(self.np)
