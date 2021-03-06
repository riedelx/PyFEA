import numpy as np
import pandas as pd
import math
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
    def __init__(self, ID, fc1, length, epsilon_t2 = 0.001, fc2_factor = 0.1, ft_factor = 1, characteristic = True,plotting=True,title="con1",Ec2 = '',Et2 = '',strain_prec=6,legend=True):
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

    def adaptic_print(self): return self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,self.ft, self.Et2

    def stmdl2_print(self):return self.Ec1,self.Ec2,-self.epsilon_2c,-self.fc2,self.Et1,self.Et2,self.epsilon_2t

    def data_frame(self,alpha=''):
        if alpha == '':
            alphac=-self.alpha
            alphat=self.alpha
        else:
            alphac=alpha[0]
            alphat=alpha[1]
        data = np.array([[self.ID, self.length, self.fc1, self.fc2, self.ft, self.Ec0, self.Ec1,
                          self.Ec2, self.Et1, self.Et2, self.Gf, self.Gc, self.epsilon_1c,
                          self.epsilon_2c, self.epsilon_1t, self.epsilon_2t, alphac, alphat]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
                      '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
                      '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
                      '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha_{c}$$', '$$alpha_{t}$$']
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
    def __init__(self,ec0,muec1,strnc1,stresc1,et0,muet1,strnt1,alphac,alphat,pseto=0,crkso=0):#pseto,crkso,
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
        self.alphac=alphac
        self.alphat=alphat

    @classmethod # alternative constructor
    def from_ADAPTIC(cls, ec1,fc1,ec2,fc2,et1,ft,et2,alphac,alphat):
        strnc1=-fc1/ec1+(fc1-fc2)/ec2
        strnt1=ft/et1-ft/et2
        return cls(ec0=ec1,muec1=ec2,strnc1=strnc1,stresc1=-fc2,et0=et1,muet1=et2,strnt1=strnt1,alphac=alphac,alphat=alphat)

    def y_E1(self,x,E1,E0,x0,x1,y0,printing=True):
        a=(E0-E1)/2/(x0-x1)
        b=E0-2*a*x0
        c=y0-a*x0**2-b*x0
        if printing: print('y(x) = {0}x**2+{1}x+{2}'.format(a,b,c))
        return (x-x0)**2*(E0-E1)/(2*(x0-x1))+E0*(x-x0)+y0

    def y_S(self,x,S,E0,x0,x1,y0,printing=True):
        a=2*(E0-S)/2/(x0-x1)
        b=E0-2*a*x0
        c=y0-a*x0**2-b*x0
        if printing: print('y(x) = {0}x**2+{1}x+{2}'.format(a,b,c))
        return (x-x0)**2*(E0-S)/(x0-x1)+E0*(x-x0)+y0

    def y_prime_E1(self,x,E1,E0,x0,x1,y0,printing=True):
        a=(x*(E0-E1)+x0*E1-x1*E0)/(x0-x1)
        b=(x**2*(E1-E0)-x0**2*(E0+E1)+2*x0*x1*E0)/(2*(x0-x1))+y0
        if printing: print('y\'(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def y_prime_S(self,x,S,E0,x0,x1,y0,printing=True):
        a=(x*2*(E0-S)+x0*(2*S-E0)-x1*E0)/(x0-x1)
        b=(S*x**2-S*x0**2-x**2*E0+x0*x1*E0)/(x0-x1)+y0
        if printing: print('y\'(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def secant(self,S,x0,y0,printing=True):
        a=S
        b=y0-a*x0
        if printing: print('secant: y(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def stress(self,strn,printing=False,retn='stress'):
        self.strn=strn

        ec0=self.ec0 # Secant compressive stiffness
        muec1=self.muec1 # Compressive softening stiffness
        strnc1=self.strnc1 # strain at residual compressive strength
        stresc1=self.stresc1 # residual compressive strength
        et0=self.et0 # Tensile stiffness
        muet1=self.muet1 # Tensile softening stiffness
        strnt1=self.strnt1 # strain at peak tensile strength ?
        pseto=self.pset # plastic strain in compression
        crkso=self.crks # plastic strain in tension
        alphac=self.alphac
        alphat=self.alphat

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain

        if(strn==pseto):
            stres=0.0
            pset=pseto
            crks=crkso
            etan=0

        if(strn<=pseto):     # Compressive strain increment

            # NOTE: alphac is the relative difference between the initial
            #       compressive tangent modulus and the secant modulus

            ec0t=(1+np.abs(alphac))*ec0
            if(alphac!=0): strnc0=(stresc1-muec1*strnc1)/(ec0-muec1)

            # Obtain the stress assuming elastic conditions, and determine
            # the force from the limiting curve

            strese=ec0t*(strn-pseto)

            if(alphac!=0 and strn>strnc0):
                stresl=strn*(ec0t+(ec0-ec0t)*strn/strnc0)
                etan=ec0t+2*(ec0-ec0t)*strn/strnc0

            elif(strn>strnc1):
                dstrn1=strn-strnc1

                # linear softening branch
                stresl=stresc1+muec1*dstrn1
                etan=muec1

                if(alphac<0):
                    # overaly with cubic function
                    dstrn0=strnc0-strnc1
                    dstrn1=dstrn1/dstrn0
                    stresl=stresl+2*alphac*muec1*dstrn0*dstrn1*(dstrn1-1)*(dstrn1-0.5)
                    etan=etan+alphac*muec1*(1-6*dstrn1*(1-dstrn1))

            else:
                # residual compressive strength
                stresl=stresc1
                etan=0.0

            # Establish the stress and the plastic set

            if(strese>stresl):
                stres=strese
                pset=pseto
                etan=ec0t

            else:
                stres=stresl
                pset=strn-stresl/ec0t

            crks=crkso

        elif(et0==0.0 or strn>=pseto+strnt1 or crkso>=strnt1):   # No tensile resistance
            stres=0.0
            pset=pseto
            crks=strnt1
            etan=0

        else:                                         # Tensile strain increment
            # NOTE: The tensile response is modified so that unloading points
            #       towards the origin of compressive loading response (i.e.
            #       plastic compressive strain), and the cracking strain is
            #       now defined as the maximum strain relative to the origin
            #       (rather than the unloading strain)

            pset=pseto

            strnt0=-muet1*strnt1/(et0-muet1)

            # Obatin relevant tensile strain to establish current stress on
            # loading/softening envelope

            dstrn=strn-pseto
            onEnv=True

            if(dstrn<=crkso): # unloading curve
                dstrn=crkso
                onEnv=False      # Use secant stiffness

            # Obtain relevant stress and tangent modulus on envelope

            if(dstrn<=strnt0):   # loading envelope
                if(alphat<=0):
                    # linear envelope: simple linear case, no need for further checks
                    stres=et0*dstrn # before et0*(strn-pseto)
                    stres=et0*(strn-pseto)
                    etan=et0
                    crks=crkso

                else:
                    # quadratic envelope
                    et0t=(1+alphat)*et0
                    stres=dstrn*(et0t+(et0-et0t)*dstrn/strnt0)
                    if(onEnv): etan=et0t+2*(et0-et0t)*dstrn/strnt0

            else:                     # softening envelope

                dstrn1=dstrn-strnt1

                # linear softening branch
                stres=muet1*dstrn1
                if(onEnv): etan=muet1

                if(alphat!=0):

                    dstrn0=strnt0-strnt1
                    dstrn1=dstrn1/dstrn0

                    if(alphat>0):
                        # overlay linear with cubic function
                        stres=stres-2*alphat*muet1*dstrn0*dstrn1*(dstrn1-1)*(dstrn1-0.5)
                        if(onEnv): etan=etan-alphat*muet1*(1-6*dstrn1*(1-dstrn1))
                    else:
                        # overlay linear with quadratic function
                        stres=stres-alphat*muet1*dstrn0*dstrn1*(dstrn1-1)
                        if(onEnv): etan=etan-alphat*muet1*(2*dstrn1-1)

            crks=dstrn

            # Determine tangent modulus as secant unloading stiffness to
            # compression origin if stress state is not on envelope

            if(not onEnv):

                etan=stres/dstrn # because dstrn=crkso
                stres=etan*(strn-pseto)

        self.strn=strn
        self.stres=stres
        self.etan=etan
        self.crks=crks
        self.pset=pset

        if retn=='etan': return self.etan
        elif retn=='all': return self.stres,self.etan
        else: return self.stres

    def plot(self,strain,retn='stress',title='con1',lineType='-',legend=True,lbl='stmdl2',xlim=(None,None),ylim=(None,None),ylabel='Stress [MPa]',xlabel='Strain',pseto='',crkso=''):
        # strain=np.arange(-np.absolute(self.strnc1),np.absolute(self.strnt1),0.0001)
        fig,ax = utils.plotBase()
        stress,etan=[],[]
        for j,i in enumerate(strain):
            if crkso!='': self.crks=crkso
            if pseto!='': self.pset=pseto
            self.stress(i)
            stress.append(self.stres)
            etan.append(self.etan)
#             if j>2:
#                 X=strain[-3:-1]
#                 Y=stress[-3:-1]
#                 print('step: {0}, etan: {1}, slope_intercept: {2}'.format(j,self.etan,self.slope_intercept(X,Y)[0]))
        stress=np.array(stress)
        etan=np.array(etan)
        strain=strain.reshape(len(strain),1)
        stress=stress.reshape(len(stress),1)
        etan=etan.reshape(len(etan),1)
        self.df=pd.DataFrame(np.hstack((stress,strain,etan)),columns=['stress','strain','etan'])
        if 'etan' in retn:
            if retn=='etan' or retn=='etan1':ax.plot(self.df['strain'],self.df['etan'],lineType, linewidth=2, markersize=5,label=lbl)
            if retn=='etan' or retn=='etan2':
                slope=self.slope_intercept(strain,stress)
                ax.plot(strain[:-1],slope,lineType, linewidth=2, markersize=5,label='slope')
        else: ax.plot(self.df['strain'],self.df['stress'],lineType, linewidth=2, markersize=5,label=lbl)
        if legend: ax.legend()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()

    def slope_intercept(self,X,Y):
        slope=[]
        for i in range(0,len(X)-1):
            if (X[i+1] - X[i])==0:
                print('ERROR: slope_intercept({},{})'.format([X[i+1],X[i]],[Y[i+1],Y[i]]))
                slope.append(math.nan)
            else: slope.append((Y[i+1] - Y[i]) / (X[i+1] - X[i]) )
        return slope

# class stmdl2_SS:
#     # this is con1 ADAPTIC model
#     def __init__(self,ec0,muec1,strnc1,stresc1,et0,muet1,strnt1,pset,crks,quad):#pseto,crkso,
#         # This subroutine calculates the stress at a monitoring point for
#         # material MODEL(2).
#
#         # Establish the stress depending on the sign of the applied
#         # strain relative to the initial plastic set and crack strain
#         self.ec0=ec0 # Secant compressive stiffness
#         self.muec1=muec1 # Compressive softening stiffness
#         self.strnc1=strnc1 # strain at residual compressive strength
#         self.stresc1=stresc1 # residual compressive strength
#         self.et0=et0 # Tensile stiffness
#         self.muet1=muet1 # Tensile softening stiffness
#         self.strnt1=strnt1 # strain at peak tensile strength ?
#         self.pseto=pseto # plastic strain in compression at the start of the step,
#                          # represents the intersection of the unloading branch with the strain axis
#         self.crkso=crkso # plastic strain in tension at the start of the step
#         self.pset=pseto
#         self.crks=crkso
#         self.quad=quad
#
#     def stress(self, strn):
#         self.strn=strn
#         if(self.strn<=self.pseto):      # Compressive strain increment
#             # NOTE: quad is the relative difference between the initial
#             #       compressive tangent modulus and the secant modulus
#
#             self.ec0t=(1+self.quad)*self.ec0 # initial tangent modulus in compression
#
#             if(self.quad>0): # implies quadratic initial compressive response
#                 self.strnc0=(self.stresc1-self.muec1*self.strnc1)/(self.ec0-self.muec1)
#
#             # Obtain the stress assuming elastic conditions, and determine
#             # the force from the limiting curve
#             self.strese=self.ec0t*(self.strn-self.pseto) # elastic stress based on ec0t
#             # this gives stress in intial compressive response if quad = 0
#
#             if(self.quad>0 and self.strn>self.strnc0): # initial quadratic response
#                 # pseto>=strn>strnc0
#
#                 #quadratic formulation for stress:
#                 self.stresl=self.strn*(self.ec0t+(self.ec0-self.ec0t)*self.strn/self.strnc0)
#                 self.etan=self.ec0t+2*(self.ec0-self.ec0t)*self.strn/self.strnc0
#
#             elif(strn>strnc1): # softening branch
#
#                 self.stresl=self.stresc1+self.muec1*(self.strn-self.strnc1)
#                 self.etan=self.muec1
#
#             else: # residual compressive strength
#
#                 self.stresl=self.stresc1
#                 self.etan=0.0
#
#             # Establish the stress and the plastic set
#
#             if(self.strese>self.stresl):
#
#                 self.stres=self.strese
#                 self.pset=self.pseto
#                 self.etan=self.ec0t
#
#             else:
#
#                 self.stres=self.stresl
#                 self.pset=self.strn-self.stresl/self.ec0t
#
#             self.crks=self.crkso
#
#         elif(self.et0==0.0 or self.strn<self.crkso+self.pseto):   # Cracked zone
#             self.stres=0.0
#             self.pset=self.pseto
#             self.crks=self.crkso
#             self.etan=0
#
#         else:   # Tensile strain increment
#
#             # Obtain the stress assuming elastic conditions, and
#             # determine the force from the limiting curve
#
#             self.strese=self.et0*(self.strn-(self.pseto+self.crkso))
#
#             self.stresl=self.muet1*(self.strn-(self.pseto+self.strnt1))+self.et0*self.strnt1 # my modification
#
#             if(self.stresl>0.0):
#
#                 self.etan=self.muet1
#
#             else:
#
#                 self.stresl=0.0
#                 self.etan=0.0
#
#             # Establish the stress and the plastic set
#
#             if(self.strese<self.stresl): # initial elastic tensile response
#
#                 self.stres=self.strese
#                 self.etan=self.et0
#                 self.crks=self.crkso
#
#             else: # tensile softening
#
#                 self.stres=self.stresl
#                 self.crks=self.strn-self.pseto-self.stresl/self.et0
#
#             self.pset=self.pseto
#         return self.stres
#
#     def plot(self,title='stmdl2',step_size=0.001,legend=True,ranges=''):
#         if ranges=='':ranges=[-np.absolute(self.strnc1),np.absolute(self.strnt1)]
#         strain=np.arange(ranges[0],ranges[1],step_size)
#         stress=np.array([self.stress(i) for i in strain])
#         strain=strain.reshape(len(strain),1)
#         stress=stress.reshape(len(stress),1)
#         self.df=pd.DataFrame(np.hstack((stress,strain)),columns=['stress','strain'])
#         self.np=plotStress(self,self.df,lbl="stmdl2",plotting='True',title=title,legend=legend)
#         self.np=np.array(self.np)
