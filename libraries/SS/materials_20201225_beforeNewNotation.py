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
        ax.set_xlabel('Strain [ε]')
        ax.set_ylabel('Stress [MPa]')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()
    return curve['strain'],curve['stress']

class con1:
    def __init__(self,ec0,muec1,strnc1,stresc1,et0,muet1,strnt1,alphac,alphat,pseto=0,crkso=0,strain_prec=6):#pseto,crkso,
        # This subroutine calculates the stress at a monitoring point for
        # material MODEL(2).

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain
        self.ec0=ec0 # Secant compressive stiffness
        self.muec1=muec1 # Compressive softening secant stiffness
        self.strnc1=strnc1 # strain at residual compressive strength
        self.stresc1=stresc1 # residual compressive strength
        self.et0=et0 # Secant tensile stiffness
        self.muet1=muet1 # Tensile softening secant stiffness
        self.strnt1=strnt1 # strain when tensile stress reaches 0 ?
        self.pset=pseto # plastic strain in compression
        self.crks=crkso # plastic strain in tension
        self.pset=pseto
        self.crks=crkso
        self.alphac=alphac
        self.alphat=alphat

        # Derived, not used in stress
        self.ec0t=round((1+np.abs(alphac))*ec0)
        self.strnc0=round((stresc1-muec1*strnc1)/(ec0-muec1),strain_prec)
        self.stresc0=round(self.strnc0*self.ec0,1)
        self.strnt0=round(-muet1*strnt1/(et0-muet1),strain_prec)
        if alphat>0: self.et0t=round((1+np.abs(alphat))*et0)
        else: self.et0t=et0
        self.ft=round(self.et0*self.strnt0,1) # Secant compressive stiffness

        data = np.array([['stmdl2',self.stresc0, self.stresc1, self.ft, self.ec0t, self.ec0,
                          self.muec1, self.et0, self.muet1, self.strnc0,
                          self.strnc1, self.strnt0, self.strnt1, self.alphac, self.alphat]])
        self.prop = pd.DataFrame(data,index=data[:,0])
        self.prop.columns = ['$$ID$$','$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$','$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$',
            '$$E_{t1}[MPa]$$','$$E_{t2}[MPa]$$','$$e_{c1}$$','$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha_{c}$$', '$$alpha_{t}$$']

    @classmethod # alternative constructor
    def from_ADAPTIC(cls, ec1,fc1,ec2,fc2,et1,ft,et2,alphac,alphat,strain_prec=6):
        strnc1=round(-fc1/ec1+(fc1-fc2)/ec2,strain_prec)
        strnt1=round(ft/et1-ft/et2,strain_prec)
        return cls(ec0=ec1,muec1=ec2,strnc1=strnc1,stresc1=-fc2,et0=et1,muet1=et2,strnt1=strnt1,alphac=alphac,alphat=alphat,strain_prec=strain_prec)

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
        pseto=self.pset # plastic strain in compression at the start of the step, represents the intersection of the unloading branch with the strain axis
        crkso=self.crks # plastic strain in tension at the start of the step, represents the intersection of the unloading branch with the strain axis
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

    def basicCurve(self,plotting=False,lbl='stmdl2',title='con1',legend=False,xlim=(None,None),ylim=(None,None),spacing=0.0001):
        strain=np.arange(-1.1*np.absolute(self.strnc1),1.1*np.absolute(self.strnt1),spacing)
        self.plot(strain,retn='stress',title=title,lineType='-',legend=legend,lbl='stmdl2',pseto=0,crkso=0,xlim=xlim,ylim=ylim)

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

class conGen(con1): # concrete material properties generator
    def __init__(self, fc1, length, ID='concrete', fc2_factor = 0.05, characteristic = True,epsilon_2t='',alpha_t='',alpha_c='',alphaT='',plotting=True,title="con1",strain_prec=6,legend=True):
        fc1 = round(fc1, 1)
        fc2 = round(fc2_factor * fc1, 1)
        if characteristic:
            fcm = round(fc1+8,2)
        else:
            fcm = fc1
        if fc1 <= 50:
            ft = round(0.3 * fcm ** (2/3), 1)
        else:
            ft = round(2.12*np.log(1+0.1*fcm), 1)
        Gf = round(73 * fcm**0.18/1000, 3)
        Ec0 = round(int(21500*(fcm/10)**(1/3)),-2)
        poisson = 0.2
        Gc = round(250 * Gf, 1)
        epsilon_1c = round(5 * fc1 / Ec0 /3, strain_prec)
        Ec1 = int(round(fc1 / epsilon_1c, -2))
        alphaC = min(max(0,round((Ec0 - Ec1)/Ec1,2)),1)
        if alpha_c=='' or alpha_c=='P' or alpha_c=='p': alphaC = alphaC
        elif alpha_c=='N' or alpha_c=='n': alphaC = -alphaC
        else: alphaC=alpha_c
        if alpha_t=='' or alpha_t=='P' or alpha_t=='p': alphaT = alphaC
        elif alpha_t=='N' or alpha_t=='n': alphaT = -alphaC
        else: alphaT=alpha_t
        epsilon_2c = round(epsilon_1c + 3 * Gc / (2 * length * fc1), strain_prec)
        Ec2 = - int(round((1 - fc2_factor) * fc1 /(epsilon_2c - epsilon_1c), -2))
        Et1 = Ec0
        epsilon_1t = round(ft / Et1, strain_prec)
        if epsilon_2t=='': epsilon_2t = round(Gf/(length*ft), strain_prec)
        else: epsilon_2t = epsilon_2t
        Et2 = - int(round(ft /(epsilon_2t - epsilon_1t), -2))

        super().__init__(Ec1,Ec2,-epsilon_2c,-fc2,Et1,Et2,epsilon_2t,alphaC,alphaT)

        data = np.array([[ID, length, fc1, fc2, ft, Ec0, Ec1,
                          Ec2, Et1, Et2, Gf, Gc, epsilon_1c,
                          epsilon_2c, epsilon_1t, epsilon_2t, alphaC,alphaT]])
        self.genProp = pd.DataFrame(data,index=data[:,0])
        self.genProp.columns = ['$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
                      '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
                      '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
                      '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$','$$alpha_{C}$$','$$alpha_{N}$$']
        self.conGen_adaptic= [Ec1, fc1, Ec2, fc2, Et1, ft, Et2]
        self.conGen_stmdl= [Ec1,Ec2,-epsilon_2c,-fc2,Et1,Et2,epsilon_2t]

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

        self.np=plotStress(self,self.df,lbl="esb1",title=title,xlim=(None,None),ylim=(None,None),plotting=plotting)
        self.np=np.array(self.np)

    def data_frame(self):
        data = np.array([[self.ID, self.fu, self.epsilon_u1, self.epsilon_u2]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$f_{u1}[MPa]$$', '$$e_{u1}$$', '$$e_{u2}$$']

class bond:
    def __init__(self,c,f_cm,ft,L,dia,n_bars,case=0,redFact=1):
        # case - refer to Table 6.1-1 MC2010:
        # 0 - Marti
        # 1 - Pull-out, good bond
        # 2 - Pull-out, all other bond cond
        # 3 - Splitting, good bond cond, unconfined
        # 4 - Splitting, good bond cond, stirrups
        # 5 - Splitting, all other bond cond, unconfined
        # 6 - Splitting, all other bond cond, stirrups
        self.case=case
        self.f_cm = f_cm # MPa
        self.ft=ft
        self.L=L
        self.dia=dia
        self.n_bars=n_bars
        self.redFact=redFact
        if self.case ==0:
            self.tau_max = self.ft * redFact
            self.s_1 = 0.001 # mm
            self.s_2 = 2 # mm
            self.s_3 = c # mm, clear distance between ribs
            self.tau_bf = self.ft * redFact
        elif self.case ==1:
            self.tau_max = 2.5 * self.f_cm**0.5 * redFact
            self.s_1 = 1 # mm
            self.s_2 = 2 # mm
            self.s_3 = c # mm, clear distance between ribs
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_max
        elif self.case ==4:
            self.tau_max = 2.5 * self.f_cm**0.5
            self.tau_bu_split=8*(f_cm/25)**0.25 * redFact
            self.s_1 = 1/self.tau_max*self.tau_bu_split # mm
            #self.s_2 = self.s_1 # mm
            self.s_3 = 0.5 * c # mm
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_bu_split
        else: print('Case error')
    def slip2tau(self,s):
        if self.case ==1:
            if 0 <= s <= self.s_1:
                tau = self.tau_max*(s/self.s_1)**self.alpha
            elif self.s_1 <= s <= self.s_2:
                tau = self.tau_max
            elif self.s_2 <= s <= self.s_3:
                tau = self.tau_max - (self.tau_max-self.tau_bf)*(s-self.s_2)/(self.s_3-self.s_2)
            else:
                tau = self.tau_bf
        elif self.case ==0:
            if 0 <= s <= self.s_1:
                tau = self.tau_max #*(s/self.s_1)
            else:
                tau = self.tau_bf
        elif self.case ==4:
            if 0 <= s <= self.s_1:
                tau = self.tau_bu_split*(s/self.s_1)**self.alpha
            elif self.s_1 <= s <= self.s_3:
                tau = self.tau_bu_split - (self.tau_bu_split-self.tau_bf)*(s-self.s_1)/(self.s_3-self.s_1)
            else:
                tau = self.tau_bf
        return tau
    def force2tau(self,force):
        return force/(self.n_bars*self.dia*np.pi*self.L)
    def tau2force(self,tau):
        return tau*(self.n_bars*self.dia*np.pi*self.L)
    def curve(self,stop=9,num=50,title='bond stress–slip relationship'):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[self.slip2tau(i) for i in x]
        ax.plot(x,y,'-', linewidth=2, markersize=5)
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Stress [MPa]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
    def astr_curve(self):
        disp1 = self.s_1
        disp2 = self.s_3
        disp3 = self.s_3+1
        f1 = self.tau2force(self.slip2tau(disp1))
        f2 = self.tau2force(self.slip2tau(disp2))
        f3 = self.tau2force(self.slip2tau(disp3))
        #np.array([[disp1, f1],[disp2, f2],[disp3, f3]]
        return f1,f2,f3,disp1,disp2,disp3
    def dataframe(self):
        return pd.DataFrame([[self.s_1,self.s_2,self.s_3,self.slip2tau(self.s_1),self.slip2tau(self.s_2),self.slip2tau(self.s_3)]],columns=['s_1','s_2','s_3','t_1','t_2','t_3'])
    def curve_force(self,stop=9,num=50,title='bond force–slip relationship'):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[(self.tau2force(self.slip2tau(i)))/1000 for i in x]
        if self.case ==0: ax.plot(x,y,'-', linewidth=2, markersize=5,label='ft')
        else: ax.plot(x,y,'-', linewidth=2, markersize=5,label='MC2010')
        f1,f2,f3,disp1,disp2,disp3=self.astr_curve()
        ax.plot([0,disp1,disp2,disp3],[0,f1/1000,f2/1000,f3/1000],'-', linewidth=2, markersize=5,label='ASTR')
        ax.legend()
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Force [kN]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
