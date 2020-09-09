import numpy as np
from numpy.linalg import inv
from numpy.linalg import multi_dot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
from scipy.integrate import quad
import sys
sys.path.insert(1, 'libraries')
import utils

# Notation
# x - positive to right
# y - positive upwards
# rz - positive anticlockwise

def k_bar(A,E,L):
    return A*E/L*np.array([[1,-1],[-1,1]])

def k_beam(I,E,L):
    #return I*E*np.array([[12/L**3,6/L**2,-12/L**3,6/L**2],[6/L**2,4/L,-6/L**2,2/L],[-12/L**3,6/L**2,12/L**3,-6/L**2],[6/L**2,2/L,-6/L**2,4/L**2]])
    return I*E/L**3*np.array([[12,6/L**2,-12,6*L],[6*L,4*L**2,-6*L,2*L**2],[-12,6*L,12,-6*L],[6*L,2*L**2,-6*L,4*L**2]])

def k_beamcol(A,I,E,L):
    return I*E/L**3*np.array([[A*L**2/I,0,0,-A*L**2/I,0,0],[0,12,6*L,0,-12,6*L],[0,6*L,4*L**2,0,-6*L,2*L**2],[-A*L**2/I,0,0,A*L**2/I,0,0],[0,-12,-6*L,0,12,-6*L],[0,6*L,2*L**2,0,-6*L,4*L**2]])

def elmCoord(elements,nodes,i):
    idx1=utils.df_index(nodes,elements['n1'][i],'name')
    idx2=utils.df_index(nodes,elements['n2'][i],'name')
    x1=nodes['x'][idx1]
    y1=nodes['y'][idx1]
    x2=nodes['x'][idx2]
    y2=nodes['y'][idx2]
    return x1,x2,y1,y2

class base2d:
    def __init__(self, materials,sections,nodes,elements,restraints,forces):
        pass

class truss2d(base2d):
    """2d truss element"""

    def __init__(self, materials,sections,nodes,elements,restraints,forces):
        super(base2d, self).__init__()
        self.materials=pd.DataFrame(materials,columns=['name','E'])
        self.sections=pd.DataFrame(sections,columns=['name','A'])
        nodes=pd.DataFrame(nodes,columns=['name','x','y'])
        self.nodes=nodes.astype({'x':'float64','y':'float64'})
        self.elements=pd.DataFrame(elements,columns=['name','n1','n2','section','material'])
        self.forces=pd.DataFrame(forces,columns=['node','x','y'])
        self.restraints=pd.DataFrame(restraints,columns=['node','x','y'])
        # add element lengths based on point coordinates
        L,x1,x2,y1,y2=[],[],[],[],[]
        for i in range(len(self.elements)):
            #if self.elements['type'][i] == 'truss':
                x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
                L.append(((x2-x1)**2+(y2-y1)**2)**0.5)
        self.elements.insert(5, "L", pd.DataFrame(L), True)

    def solve(self):
        # Assemble element and global stiffness matrices
        k_glob = np.zeros((len(self.nodes)*2,len(self.nodes)*2))
        k_elms =[]
        G_elms=[]
        indices=[]
        for i in range(len(self.elements)):
            #if self.elements['type'][i] == 'truss':
                A = utils.df_value(self.sections,self.elements['section'][i],'name','A')
                E = utils.df_value(self.materials,self.elements['material'][i],'name','E')
                idx1=utils.df_index(self.nodes,self.elements['n1'][i],'name')
                idx2=utils.df_index(self.nodes,self.elements['n2'][i],'name')
                L=self.elements['L'][i]
                x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
                s=(y2-y1)/self.elements['L'][i]
                c=(x2-x1)/self.elements['L'][i]
                G=np.array([[c,s,0,0],[0,0,c,s]])
                G_elms.append(G)
                k_elm=multi_dot([np.transpose(G),k_bar(A,E,L),G])
                k_elms.append(k_elm)
                k_glob=k_glob.reshape(k_glob.shape[0],k_glob.shape[1])
                k_elm=k_elm.reshape(k_elm.shape[0],k_elm.shape[1])
                indices.append([idx1*2,idx1*2+1,idx2*2,idx2*2+1]) # corresponding indices of the element nodes
                k_glob[np.ix_(indices[i],indices[i])]=k_glob[np.ix_(indices[i],indices[i])]+k_elm

        # Assemble external force vector
        f_ext = np.zeros((len(self.nodes)*2,1))
        for i in range(len(self.forces)):
            #if self.elements['type'][i] == 'truss':
                idx=utils.df_index(self.nodes,self.forces['node'][i],'name')
                f_ext[idx*2]=f_ext[idx*2]+self.forces['x'][i]
                f_ext[idx*2+1]=f_ext[idx*2+1]+self.forces['y'][i]
        # free DoF
        restdof = []
        for i in range(len(self.restraints)):
            idx=utils.df_index(self.nodes,self.restraints['node'][i],'name')
            if self.restraints['x'][i]==1:
                restdof.append(idx*2)
            if self.restraints['y'][i]==1:
                restdof.append(idx*2+1)
        freedof=np.arange(len(self.nodes)*2)
        freedof=np.delete(freedof,restdof)
        # Solve displacements
        self.disp = np.zeros((len(self.nodes)*2,1))
        self.disp[freedof]=multi_dot([inv(k_glob[np.ix_(freedof,freedof)]),f_ext[freedof]])
        # Element force vector (external and reaction forces on nodes)
        f_efv = multi_dot([k_glob,self.disp])
        # Reaction forces
        self.f_reac = f_efv - f_ext

        # Post-processing
        # displaced coordinates of the nodes
        # displaced_nodes=self.nodes.copy()
        # for i in range(len(self.nodes)):
        #     displaced_nodes.at[i,'x']=self.nodes['x'][i]+self.disp[i*2]
        #     displaced_nodes.at[i,'y']=self.nodes['y'][i]+self.disp[i*2+1]

        # axial forces
        f_elm=[]
        f_int = []
        for i in range(len(self.elements)):
            f=multi_dot([k_elms[i],self.disp[indices[i]]])
            f_int.append(f)
            f=multi_dot([G_elms[i],f])
            f_elm.append(f)
        f_elm=np.array(f_elm).flatten()
        self.axial=[f_elm[i*2] for i in range(int(len(f_elm)/2))]

    def plot_deformed(self,disp_magn=10,title='deformed shape',markersize=10,grid=True,figsize=(6,4)):
        displacements=self.nodes.copy()
        for i in range(len(self.nodes)):
            displacements.at[i,'x']=self.nodes['x'][i]+self.disp[i*2]*disp_magn
            displacements.at[i,'y']=self.nodes['y'][i]+self.disp[i*2+1]*disp_magn
        fig,ax = utils.plotBase(grid=grid,figsize=figsize)
        for i in range(len(self.elements)):
            x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
            ax.plot([x1,x2],[y1,y2],'-', linewidth=2, markersize=5,color='b')
        for i in range(len(self.elements)):
            x1,x2,y1,y2 = elmCoord(self.elements,displacements,i)
            ax.plot([x1,x2],[y1,y2],'-', linewidth=2, markersize=5,color='r')
        for i in range(len(self.restraints)):
            x=utils.df_value(self.nodes,self.restraints['node'][i],'name','x')
            y=utils.df_value(self.nodes,self.restraints['node'][i],'name','y')
            if self.restraints['x'][i]==1 and self.restraints['y'][i]==1:
                ax.plot(x,y,'s', linewidth=2, markersize=markersize,color='k')
            if self.restraints['x'][i]==1 and self.restraints['y'][i]==0:
                ax.plot(x,y,'>', linewidth=2, markersize=markersize,color='k')
            if self.restraints['x'][i]==0 and self.restraints['y'][i]==1:
                ax.plot(x,y,'^', linewidth=2, markersize=markersize,color='k')
        ax.set_title(title)
        plt.axis('equal')
        plt.show()

    def plot_var(self,var,disp_magn=10,title='variable',markersize=10,grid=True,scale=1,figsize=(6,4)):
        fig,ax = utils.plotBase(grid=grid,figsize=figsize)
        var = np.array(getattr(self, var))*scale
        norm = plt.Normalize(np.min(var), np.max(var))
        cmap = plt.get_cmap('seismic')
        c = cmap(norm(var))
        for i in range(len(self.elements)):
            x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
            ax.plot([x1,x2],[y1,y2],'-', linewidth=2, markersize=5,c=c[i])
        for i in range(len(self.restraints)):
            x=utils.df_value(self.nodes,self.restraints['node'][i],'name','x')
            y=utils.df_value(self.nodes,self.restraints['node'][i],'name','y')
            if self.restraints['x'][i]==1 and self.restraints['y'][i]==1:
                ax.plot(x,y,'s', linewidth=2, markersize=markersize,color='k')
            if self.restraints['x'][i]==1 and self.restraints['y'][i]==0:
                ax.plot(x,y,'>', linewidth=2, markersize=markersize,color='k')
            if self.restraints['x'][i]==0 and self.restraints['y'][i]==1:
                ax.plot(x,y,'^', linewidth=2, markersize=markersize,color='k')
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, ticks=var)
        plt.axis('equal')
        ax.set_title(title)
        plt.show()

class beamcol2d:
    """2d beam-column element"""

    def __init__(self, materials,sections,nodes,elements,restraints,forces_nodes,forces_distributed):
        #super(, self).__init__()
        self.materials=pd.DataFrame(materials,columns=['name','E'])
        self.sections=pd.DataFrame(sections,columns=['name','A','I'])
        nodes=pd.DataFrame(nodes,columns=['name','x','y'])
        self.nodes=nodes.astype({'x':'float64','y':'float64'})
        self.elements=pd.DataFrame(elements,columns=['name','n1','n2','section','material'])
        if forces_nodes != None:
            self.forces_nodes=pd.DataFrame(forces_nodes,columns=['node','x','y','rz'])
        else:
            self.forces_nodes=None
        if forces_distributed != None:
            self.forces_distributed=pd.DataFrame(forces_distributed,columns=['element','n1','n2'])
        else:
            self.forces_distributed = None
        self.restraints=pd.DataFrame(restraints,columns=['node','x','y','rz'])
        # add element lengths based on point coordinates
        L,x1,x2,y1,y2=[],[],[],[],[]
        for i in range(len(self.elements)):
            x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
            L.append(((x2-x1)**2+(y2-y1)**2)**0.5)
        self.elements.insert(5, "L", pd.DataFrame(L), True)

        # Assemble element and global stiffness matrices
        k_glob = np.zeros((len(self.nodes)*3,len(self.nodes)*3))
        k_elms =[]
        G_elms=[]
        indices=[]
        for i in range(len(self.elements)):
            A = utils.df_value(self.sections,self.elements['section'][i],'name','A')
            I = utils.df_value(self.sections,self.elements['section'][i],'name','I')
            E = utils.df_value(self.materials,self.elements['material'][i],'name','E')
            idx1=utils.df_index(self.nodes,self.elements['n1'][i],'name')
            idx2=utils.df_index(self.nodes,self.elements['n2'][i],'name')
            L=self.elements['L'][i]
            x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
            s=(y2-y1)/self.elements['L'][i]
            c=(x2-x1)/self.elements['L'][i]
            G=np.array([[c,s,0,0,0,0],[-s,c,0,0,0,0],[0,0,1,0,0,0],[0,0,0,c,s,0],[0,0,0,-s,c,0],[0,0,0,0,0,1]])
            G_elms.append(G)
            k_elm=multi_dot([np.transpose(G),k_beamcol(A,I,E,L),G])
            k_elms.append(k_elm)
            k_glob=k_glob.reshape(k_glob.shape[0],k_glob.shape[1])
            k_elm=k_elm.reshape(k_elm.shape[0],k_elm.shape[1])
            indices.append([idx1*3,idx1*3+1,idx1*3+2,idx2*3,idx2*3+1,idx2*3+2]) # corresponding indices of the element nodes
            k_glob[np.ix_(indices[i],indices[i])]=k_glob[np.ix_(indices[i],indices[i])]+k_elm
        self.k_glob=k_glob

        # Assemble external force vector based on nodal loads
        self.f_ext_nodes = np.zeros((len(self.nodes)*3,1))
        if forces_nodes != None:
            for i in range(len(self.forces_nodes)):
                idx=utils.df_index(self.nodes,self.forces_nodes['node'][i],'name')
                self.f_ext_nodes[idx*3]=self.f_ext_nodes[idx*3]+self.forces_nodes['x'][i]
                self.f_ext_nodes[idx*3+1]=self.f_ext_nodes[idx*3+1]+self.forces_nodes['y'][i]
                self.f_ext_nodes[idx*3+2]=self.f_ext_nodes[idx*3+2]+self.forces_nodes['rz'][i]

        # convert distributed element load to nodal loads
        self.f_ext_distr = np.zeros((len(self.nodes)*3,1))
        self.elm_f_distr=np.zeros((len(self.elements)*6,1)) # element distributed forces converted to equivalent axial, shear and moments
        if forces_distributed != None:
            f_temp=[]
            for i in range(len(self.forces_distributed)):
                idx_elm=utils.df_index(self.elements,self.forces_distributed['element'][i],'name')
                idx_n1=utils.df_index(self.nodes,self.elements['n1'][idx_elm],'name')
                idx_n2=utils.df_index(self.nodes,self.elements['n2'][idx_elm],'name')
                L=self.elements['L'][idx_elm]
                w1=self.forces_distributed['n1'][i]
                w2=self.forces_distributed['n2'][i]
                area=np.trapz([w1,w2], x=[0,L])
                centroid = utils.centroidX([[0,w1],[L,w2]])
                m1=-L**2/10*(w1/2+w2/3)
                m2=L**2/10*(w1/3+w2/2)
                r1=-(m1+m2+area*(L-centroid))/L
                r2=-area-r1
                f_temp.append([m1,m2,r1,r2])
                temp=multi_dot([G_elms[idx_elm],np.array([0,r1,m1,0,r2,m2])]).reshape(6,1)
                self.elm_f_distr[idx_elm*6:idx_elm*6+6]=self.elm_f_distr[idx_elm*6:idx_elm*6+6]+temp
                self.f_ext_distr[idx_n1*3]=self.f_ext_distr[idx_n1*3]-temp[0]
                self.f_ext_distr[idx_n1*3+1]=self.f_ext_distr[idx_n1*3+1]-temp[1]
                self.f_ext_distr[idx_n1*3+2]=self.f_ext_distr[idx_n1*3+2]-temp[2]
                self.f_ext_distr[idx_n2*3]=self.f_ext_distr[idx_n2*3]-temp[3]
                self.f_ext_distr[idx_n2*3+1]=self.f_ext_distr[idx_n2*3+1]-temp[4]
                self.f_ext_distr[idx_n2*3+2]=self.f_ext_distr[idx_n2*3+2]-temp[5]
            self.forces_distributed=self.forces_distributed.join(pd.DataFrame(f_temp,columns=['m1','m2','r1','r2']))

        # comnine nodal and distributed loads
        self.f_ext = self.f_ext_nodes + self.f_ext_distr

        # free DoF
        restdof = []
        for i in range(len(self.restraints)):
            idx=utils.df_index(self.nodes,self.restraints['node'][i],'name')
            if self.restraints['x'][i]==1:
                restdof.append(idx*3)
            if self.restraints['y'][i]==1:
                restdof.append(idx*3+1)
            if self.restraints['rz'][i]==1:
                restdof.append(idx*3+2)
        freedof=np.arange(len(self.nodes)*3)
        freedof=np.delete(freedof,restdof)
        # Solve displacements
        self.disp = np.zeros((len(self.nodes)*3,1))
        self.disp[freedof]=multi_dot([inv(k_glob[np.ix_(freedof,freedof)]),self.f_ext[freedof]])
        # Element force vector (external and reaction forces on nodes)
        f_efv = multi_dot([k_glob,self.disp])
        # Reaction forces
        self.f_reac = f_efv - self.f_ext

        # Post-processing
        # displaced coordinates of the nodes
        # displaced_nodes=self.nodes.copy()
        # for i in range(len(self.nodes)):
        #     displaced_nodes.at[i,'x']=self.nodes['x'][i]+self.disp[i*2]
        #     displaced_nodes.at[i,'y']=self.nodes['y'][i]+self.disp[i*2+1]

        # axial forces
        f_elm=[]
        f_int = []
        for i in range(len(self.elements)):
            f=multi_dot([k_elms[i],self.disp[indices[i]]])
            f_int.append(f)
            f=multi_dot([G_elms[i],f])
            f_elm.append(f)
        self.f_elm=np.array(f_elm).flatten()+self.elm_f_distr.reshape(len(self.elm_f_distr))
        self.axial=[self.f_elm[i*3] for i in range(int(len(self.f_elm)/3))]
        self.shear=[self.f_elm[i*3+1] for i in range(int(len(self.f_elm)/3))]
        self.moment=[self.f_elm[i*3+2] for i in range(int(len(self.f_elm)/3))]

    def plot_moments(self,title='moment diagram',figsize=(10,4),scale=1,steps=10,fontsize=12,round=3):
        scale = scale * np.max(self.elements['L'])/3 / np.max(np.abs(self.moment))
        relativeOffset = scale * np.max(self.elements['L'])
        fig, ax = plt.subplots(figsize=figsize)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        element_summary=[]
        for i in range(len(self.elements)):
            #elements
            x_el=[utils.df_value(self.nodes,self.elements['n1'][i],'name','x'),utils.df_value(self.nodes,self.elements['n2'][i],'name','x')]
            y_el=[utils.df_value(self.nodes,self.elements['n1'][i],'name','y'),utils.df_value(self.nodes,self.elements['n2'][i],'name','y')]
            ax.plot(x_el,y_el,'-', linewidth=3, marker='o',color = 'black')
            #moments
            x1,x2,y1,y2 = elmCoord(self.elements,self.nodes,i)
            s=(y2-y1)/self.elements['L'][i]
            c=(x2-x1)/self.elements['L'][i]
            if self.elements['name'][i] in self.forces_distributed['element'].values:
                x_mom=[x_el[0]]
                y_mom=[y_el[0]]
                moments=[]
                L=self.elements['L'][i]
                steps = steps
                w1 = utils.df_value(self.forces_distributed,self.elements['name'][i],'element','n1')
                w2 = utils.df_value(self.forces_distributed,self.elements['name'][i],'element','n2')
                for j in range(0,steps+1):
                    x=L/steps * j
                    w2_temp=w1+(w2-w1)*x/L
                    # area=np.trapz([w1,w2_temp], x=[0,L])
                    # leverarm = w2_temp-utils.centroidX([[0,w1],[x,w2_temp]])
                    area=(w1+w2_temp)/2*x
                    if area!=0: leverarm=(w1*x**2/2+(w2_temp-w1)*x**2/6)/area
                    else: leverarm=0
                    moment=-(-self.moment[i*2]+self.shear[i*2]*x+area*leverarm) #upsidedown moments
                    moments.append(moment)
                    x_mom.append(utils.df_value(self.nodes,self.elements['n1'][i],'name','x')+x+moment*s*scale)
                    y_mom.append(utils.df_value(self.nodes,self.elements['n1'][i],'name','y')+moment*c*scale)
                # x_mom.append(utils.df_value(self.nodes,self.elements['n2'][i],'name','x')+s*self.moment[i*2+1]*scale)
                # y_mom.append(utils.df_value(self.nodes,self.elements['n2'][i],'name','y')+c*self.moment[i*2+1]*scale)
                x_mom.append(x_el[1])
                y_mom.append(y_el[1])
                ax.text(x_mom[1]+relativeOffset*c,y_mom[1]+utils.matchSign(relativeOffset,y_mom[1]), utils.round_sig(moments[0],round), fontsize=fontsize)
                ax.text(x_mom[-2]-4*relativeOffset*c,y_mom[-2]+utils.matchSign(relativeOffset,y_mom[-2]), utils.round_sig(moments[-1],round), fontsize=fontsize)
                if utils.maxAbs(moments[0],moments[-1])>0:indice = min((val, idx) for (idx, val) in enumerate(moments))[1]
                else:indice = max((val, idx) for (idx, val) in enumerate(moments))[1]
                ax.text(x_mom[indice+1],y_mom[indice]+utils.matchSign(relativeOffset,y_mom[indice+1]), utils.round_sig(moments[indice],round), fontsize=fontsize)
                element_summary.append([self.elements['name'][i],utils.round_sig(moments[0],round),utils.round_sig(moments[-1],round),utils.round_sig(moments[indice],round),self.shear[i*2],self.shear[i*2+1]])
            else:
                x_mom=[x_el[0],utils.df_value(self.nodes,self.elements['n1'][i],'name','x')+s*self.moment[i*2]*scale,utils.df_value(self.nodes,self.elements['n2'][i],'name','x')-s*self.moment[i*2+1]*scale,x_el[1]]
                y_mom=[y_el[0],utils.df_value(self.nodes,self.elements['n1'][i],'name','y')+c*self.moment[i*2]*scale,utils.df_value(self.nodes,self.elements['n2'][i],'name','y')-c*self.moment[i*2+1]*scale,y_el[1]]
                ax.text(x_mom[1]+relativeOffset*c,y_mom[1]+utils.matchSign(relativeOffset,y_mom[1]), utils.round_sig(self.moment[i*2],round), fontsize=fontsize)
                ax.text(x_mom[-2]-4*relativeOffset*c,y_mom[-2]+utils.matchSign(relativeOffset,y_mom[-2]), utils.round_sig(self.moment[i*2+1],round), fontsize=fontsize)
                element_summary.append([self.elements['name'][i],utils.round_sig(self.moment[i*2],round),utils.round_sig(self.moment[i*2+1],round),(utils.round_sig(self.moment[i*2],round)+utils.round_sig(self.moment[i*2+1],round))/2,self.shear[i*2],self.shear[i*2+1]])
            ax.plot(x_mom,y_mom,'-', linewidth=2, marker='',color = 'blue')
        self.element_summary=pd.DataFrame(element_summary,columns=['name','m1','m2','mMid','S1','S2'])
        plt.axis('equal')
        ax.set_title(title)
        plt.show()

def textOffset(val,relativeOffset,cosSin,end=1):
    if val >= 0 and end==1: return val+relativeOffset*cossin
