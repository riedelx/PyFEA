import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import utils

def plotStress(curve,title="",xlim=(None,None),ylim=(None,None)):
    fig,ax = utils.plotBase()
    ax.plot(curve[0],curve[1],'-', linewidth=2, markersize=5)
    ax.set_title(title)
    ax.set_xlabel('Strain [ε]')
    ax.set_ylabel('Stress [MPa]')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.show()

class con1:
    def __init__(self,Ec1,fc1,Ec2,fc2,Et1,ft,Et2,alphac,alphat,pset=0,crks=0):#pset,crks,
        # This subroutine calculates the stress at a monitoring point for
        # material MODEL(2).

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain
        self.Ec1=Ec1 # secant compressive stiffness
        self.fc1=fc1 # peak compressive strength
        self.Ec2=Ec2 # compressive softening secant stiffness
        self.fc2=fc2 # residual compressive strength
        self.Et1=Et1 # secant tensile stiffness
        self.ft=ft # peak tensile strength
        self.Et2=Et2 # tensile softening secant stiffness
        self.pset=pset # plastic strain in compression
        self.crks=crks # cracking strain in tension
        self.alphac=alphac
        self.alphat=alphat
        self.adaptic = [Ec1,-fc1,Ec2,-fc2,Et1,ft,Et2,alphac,alphat]
        self.Ec0=(1+np.abs(alphac))*Ec1 # secant compressive stiffness
        self.eps_c1=fc1/Ec1 # strain at peak compressive strength
        self.eps_c2=self.eps_c1-fc1/Ec2 # strain at residual compressive strength
        self.eps_t1=ft/Et1 # strain at peak tensile strength
        self.eps_t2=self.eps_t1-ft/Et2 # strain when tensile stress reaches 0

        # # Establish the stress depending on the sign of the applied
        # # strain relative to the initial plastic set and crack strain
        # self.Ec1=Ec1 # secant compressive stiffness
        # self.Ec2=Ec2 # compressive softening secant stiffness
        # self.eps_c2=eps_c2 # strain at residual compressive strength
        # self.fc2=fc2 # residual compressive strength
        # self.Et1=Et1 # secant tensile stiffness
        # self.Et2=Et2 # tensile softening secant stiffness
        # self.eps_t2=eps_t2 # strain when tensile stress reaches 0
        # self.pset=pset # plastic strain in compression
        # self.crks=crks # cracking strain in tension
        # self.alphac=alphac
        # self.alphat=alphat
        # self.input = [Ec1,Ec2,eps_c2,fc2,Et1,Et2,eps_t2,alphac,alphat]
        # # Derived, not used in stress
        # self.Ec0=(1+np.abs(alphac))*Ec1 # secant compressive stiffness
        # self.eps_c1=(fc2-Ec2*eps_c2)/(Ec1-Ec2) # strain at peak compressive strength
        # self.fc1=self.eps_c1*self.Ec1 # peak compressive strength
        # self.eps_t1=-Et2*eps_t2/(Et1-Et2) # strain at peak tensile strength
        # if alphat>0: self.et0t=(1+np.abs(alphat))*Et1 # secant tensile stiffness
        # else: self.et0t=Et1
        # self.ft=self.Et1*self.eps_t1 # peak tensile strength

        data = np.array([[type(self).__name__,self.fc1, self.fc2, self.ft, self.Ec0, self.Ec1,
                          self.Ec2, self.Et1, self.Et2, self.eps_c1,
                          self.eps_c2, self.eps_t1, self.eps_t2, self.alphac, self.alphat]])
        self.prop = pd.DataFrame(data,index=data[:,0])
        self.prop.columns = ['type','$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$','$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$',
            '$$E_{t1}[MPa]$$','$$E_{t2}[MPa]$$','$$e_{c1}$$','$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha_{c}$$', '$$alpha_{t}$$']

    # @classmethod # alternative constructor
    # def from_ADAPTIC(cls, ec1,fc1,ec2,fc2,et1,ft,et2,alphac,alphat):
    #     eps_c2=-fc1/ec1+(fc1-fc2)/ec2
    #     eps_t2=ft/et1-ft/et2
    #     return cls(Ec1=ec1,Ec2=ec2,eps_c2=eps_c2,fc2=-fc2,Et1=et1,Et2=et2,eps_t2=eps_t2,alphac=alphac,alphat=alphat)

    def stress(self,strn):
        self.strn=strn

        ec0=self.Ec1 # secant compressive stiffness
        muec1=self.Ec2 # compressive softening secant stiffness
        strnc1=self.eps_c2 # strain at residual compressive strength
        stresc1=self.fc2 # residual compressive strength
        et0=self.Et1 # secant tensile stiffness
        muet1=self.Et2 # tensile softening secant stiffness
        strnt1=self.eps_t2 # strain when tensile stress reaches 0
        pseto=self.pset # plastic strain in compression at the start of the step, represents the intersection of the unloading branch with the strain axis
        crkso=self.crks # cracking strain in tension at the start of the step, represents the intersection of the unloading branch with the strain axis
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

        # self.strn=strn
        # self.stres=stres
        # self.etan=etan
        # self.crks=crks
        # self.pset=pset

        return stres

    def basicCurve(self,title='',xlim=(None,None),ylim=(None,None),spacing=0.00001):
        strain=np.arange(-1.1*np.absolute(self.eps_c2),1.1*self.eps_t2,spacing)
        self.plot(strain,title=title,lineType='-',xlim=xlim,ylim=ylim) # pset=0,crks=0

    def plot(self,strain,title='',lineType='-',xlim=(None,None),ylim=(None,None)): #,pset=0,crks=0
        stress=[]
        for j,i in enumerate(strain):
            # if crks!='': self.crks=crks
            # if pset!='': self.pset=pset
            stress.append(self.stress(i))
        plotStress([strain,stress],title=title,xlim=xlim,ylim=ylim)

class con1gen(con1): # concrete material properties generator
    def __init__(self, fc1, length, fc2_factor = 0.05, characteristic = True, epsilon_2t='', alphac=''):
        fc2 = fc2_factor * fc1
        if characteristic: fcm = fc1+8
        else: fcm = fc1
        if fc1 <= 50: ft = 0.3 * fcm ** (2/3)
        else: ft = 2.12*np.log(1+0.1*fcm)
        Gf = 73 * fcm**0.18/1000
        Ec0 = 21500*(fcm/10)**(1/3) # initial compressive stiffness
        poisson = 0.2
        Gc = 250 * Gf
        epsilon_1c =5 * fc1 / Ec0 /3
        epsilon_2c = epsilon_1c + 3 * Gc / (2 * length * fc1)
        Ec1 = fc1 / epsilon_1c # secant compressive stiffness
        Ec2 = -(fc1-fc2)/(epsilon_2c - epsilon_1c) # secant compressive softening stiffness
        alpha = min(max(0,(Ec0 - Ec1)/Ec1),1)
        if alphac=='' or alphac=='N' or alphac=='n': alphac = -alpha
        elif alphac=='P' or alphac=='p': alphac = alpha
        else: alphac=alphac
        alphat = -1
        Ec0 = (1+np.abs(alphac))*Ec1 # in case alphac is changed
        Et1 = Ec0 # secant tensile stiffness
        epsilon_1t = ft / Et1
        # epsilon_2t
        area_f = Gf/length # fracture energy area
        area_f_soft = area_f - epsilon_1t*ft/2 # area under the softening curve
        # Figure 1 page p 17 in RTD 2010
        eps_u = 2*area_f_soft/ft
        E0 = -ft/(eps_u-epsilon_1t) # tangent stiffness at epsilon_1t
        E1 = 0 # tangent stiffness at epsilon_2t
        epsilon_2t = max((epsilon_1t*E0+epsilon_1t*E1-2*ft)/(E0+E1),epsilon_1t)
        # print('epsilon_1t={},eps_u={},epsilon_2t={}'.format(epsilon_1t,eps_u,epsilon_2t))
        # print('area_f={},area_f_soft={}'.format(area_f,area_f_soft))
        # print('E0={},E1={}'.format(E0,E1))
        # epsilon_2t using area under parabola curve
        # epsilon_2t = 3*area_t_par/ft
        Et2 = - ft /(epsilon_2t - epsilon_1t) # secant tensile softening stiffness

        super().__init__(Ec1,-fc1,Ec2,-fc2,Et1,ft,Et2,alphac,alphat)
        #super().__init__(Ec1,Ec2,-epsilon_2c,-fc2,Et1,Et2,epsilon_2t,alphac,alphat)

        # data = np.array([[type(self).__name__,length, fc1, fc2, ft, Ec0, Ec1,
        #                   Ec2, Et1, Et2, Gf, Gc, epsilon_1c,
        #                   epsilon_2c, epsilon_1t, epsilon_2t, alphac,alphat]])
        # self.genProp = pd.DataFrame(data,index=data[:,0])
        # self.genProp.columns = ['type','$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
        #               '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
        #               '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
        #               '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$','$$alpha_{C}$$','$$alpha_{T}$$']

class epm1: # elsto-platic model
    def __init__(self, E1, fy, fu, epsilon_u,tension=True,compression=True):
        self.E1 = E1
        self.fy = fy
        self.fu = fu
        self.epsilon_u = epsilon_u
        self.epsilon_y = fy / E1
        self.E2 = (fu - fy) / (epsilon_u - self.epsilon_y)
        self.mu = self.E2 / E1
        self.tension=tension
        self.compression=compression

        data = np.array([[type(self).__name__, self.E1, self.E2, self.fy, self.fu, self.epsilon_y, self.epsilon_u,
                          self.mu]])
        self.prop = pd.DataFrame(data,index=data[:,0])
        self.prop.columns = ['type', '$$E_{1}[MPa]$$', '$$E_{2}[MPa]$$', '$$f_{y}[MPa]$$', '$$f_{u}[MPa]$$',
                      '$$e_{y}$$','$$e_{u}$$','$$mu$$']

    def stress(self,strain):
        if self.tension and 0<strain<=self.epsilon_y: return strain*self.E1
        elif self.tension and self.epsilon_y<strain<=self.epsilon_u: return self.fy+(strain-self.epsilon_y)*self.E2
        elif self.compression and 0>strain>=-self.epsilon_y: return strain*self.E1
        elif self.compression and -self.epsilon_y>strain>=-self.epsilon_u: return -self.fy+(strain-self.epsilon_y)*self.E2
        else: return 0

    def basicCurve(self,title='',xlim=(None,None),ylim=(None,None)):
        # data = [[0,self.epsilon_y,self.epsilon_u],[0,self.fy,self.fu]]
        # if self.tension and self.compression: data=[[[-i for i in reversed(z)]+z] for z in data]
        # elif self.compression:data=[[[-i for i in reversed(z)]] for z in data]
        data = [0,self.epsilon_y,self.epsilon_u]
        if self.tension and self.compression: data=[-i for i in reversed(data)]+data
        elif self.compression:data=[-i for i in reversed(data)]
        self.plot(data,title=title,xlim=xlim,ylim=ylim)

    def plot(self,strain,title='',lineType='-',xlim=(None,None),ylim=(None,None)):
        stress=[self.stress(i) for i in strain]
        plotStress([strain,stress],title=title,xlim=xlim,ylim=ylim)

class esb1: # equivalent stress block
    def __init__(self, fu, lamb, epsilon_u):
        self.fu = fu
        self.epsilon_u = epsilon_u
        self.epsilon_y = lamb*self.epsilon_u
        self.lamb=lamb

        data = np.array([[type(self).__name__, self.fu, self.epsilon_u,self.lamb]])
        self.prop = pd.DataFrame(data,index=data[:,0])
        self.prop.columns = ['type', '$$f_{u1}[MPa]$$', '$$e_{u}$$', '$$\lambda$$']

    def stress(self,strain):
        if self.epsilon_y <= strain <= self.epsilon_u: return self.fu
        else: return 0

    def basicCurve(self,title='',xlim=(None,None),ylim=(None,None)):
        strain = [0,0.999*self.epsilon_y,self.epsilon_y,self.epsilon_u]
        self.plot(strain,title=title,xlim=xlim,ylim=ylim)

    def plot(self,strain,title='',lineType='-',xlim=(None,None),ylim=(None,None)):
        stress=[self.stress(i) for i in strain]
        plotStress([strain,stress],title=title,xlim=xlim,ylim=ylim)
