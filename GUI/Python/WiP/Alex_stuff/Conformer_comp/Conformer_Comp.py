# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:58:36 2024

@author: AHaack
"""

import os
import sys
import numpy as np
import numpy.linalg as LA
import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster import hierarchy

Eh2eV = 27.21139
Eh2kJ = 2625.5
subfol = 'extract'

a2m = {'H' :  1.,                                                              'He':  4.,
       'Li':  7., 'Be':  9., 'B': 11, 'C': 12., 'N': 14., 'O': 16., 'F' : 19., 'Ne': 20.,
       'Na': 23., 'Mg': 24.,                    'P': 31., 'S': 32., 'Cl': 35., 'Ar': 40.,
       } # from atom label to mass

at_cols = {'C': 'black', 'H': 'dimgrey', 'O': 'orangered', 'N': 'royalblue',
           'S': 'gold', 'F': 'cyan', 'Cl': 'lime', 'Br': 'firebrick', 'P': 'orange'}

#%% make ORCA inp files
method = ['! wB97X-D3BJ Opt Freq def2-TZVPP def2/J RIJCOSX ExtremeSCF',
          '! PAL8 VeryTightOpt defgrid3',
          '',
          '%maxcore 4096',
          '',
          '%geom',
          '  MaxIter 400',
          '  # Step qn',
          '  # MaxStep 0.1',
          '  Calc_Hess true',
          '  # Recalc_Hess 2',
          'end',
          '',
          '%method',
          '  Z_tol 1e-12',
          'end',
          '',
          '%chelpg',
          '  grid 0.1',
          '  rmax 3.0',
          'end',
          '',
          '%elprop',
          '  Dipole true',
          '  Quadrupole true',
          '  Polar 1',
          'end',
          '',
          '* xyz 1 1']

#%% functions

def rot_x(theta):
    '''return a rotation matrix around x-axis by angle theta (in rad)'''
    Rx = np.array([[1., 0.            , 0.            ],
                   [0., +np.cos(theta), -np.sin(theta)],
                   [0., +np.sin(theta), +np.cos(theta)]])
    return Rx

def rot_y(theta):
    '''return a rotation matrix around y-axis by angle theta (in rad)'''
    Ry = np.array([[+np.cos(theta), 0., +np.sin(theta)],
                   [ 0.           , 1.,     0.        ],
                   [-np.sin(theta), 0., +np.cos(theta)]])
    return Ry
    
def rot_z(theta):
    '''return a rotation matrix around z-axis by angle theta (in rad)'''
    Rz = np.array([[+np.cos(theta), -np.sin(theta), 0.],
                   [+np.sin(theta), +np.cos(theta), 0.],
                   [    0.        ,     0.        , 1.]])
    return Rz

def orient(xyz_in,fix_ats):
    '''Orientes the molecule such that the fix_ats have the following properties:
    1st one lies at (0,0,0), 2nd one lies at (+x,0,0), 3rd one lies at (x,+y,0)'''
    xyz = np.copy(xyz_in)
    ## translate pep[0] to origin
    xyz -= xyz[fix_ats[0]]
    
    ## rotate pep[1]-pep[0] onto x-axis
    # first, rotate around z such that x-component of pep[1] is positive
    if xyz[fix_ats[1],0] < 0.:
        Rz_180 = rot_z(np.pi)
        xyz = np.array([np.dot(Rz_180,xyz_i) for xyz_i in xyz])
    
    # get pep[1]-pep[0] vector
    bond = xyz[fix_ats[1]]-xyz[fix_ats[0]]
    # get rotation matrix around z such that y=0
    bond_xy = np.array([bond[0],bond[1],0.]) # projection onto xy plane
    theta_z = np.arccos(bond_xy[0]/LA.norm(bond_xy))
    theta_z *= -np.sign(bond[1]) # if sign(y)=+ -> -; else +
    Rz = rot_z(theta_z)
    # get rotation matrix around y such that z=0
    theta_y = np.abs(np.arccos(np.dot(bond,bond_xy)/(LA.norm(bond)*LA.norm(bond_xy))))
    theta_y *= np.sign(bond[2]) # if sign(z)=+ -> +; else -
    Ry = rot_y(theta_y)
    # do rotation
    xyz = np.array([np.dot(Ry,np.dot(Rz,xyz_i)) for xyz_i in xyz])
    
    ## rotate pep[2] into xy plane
    # y-component of pep[2] should be positive
    if xyz[fix_ats[2],1] < 0.:
        Rx_180 = rot_x(np.pi)
        xyz = np.array([np.dot(Rx_180,xyz_i) for xyz_i in xyz])
    # get new bond
    bond = xyz[fix_ats[2]]-xyz[fix_ats[0]]
    bond_yz = np.array([0.,bond[1],bond[2]]) # projection onto yz plane
    theta_x = np.abs(np.arccos(bond_yz[1]/LA.norm(bond_yz)))
    theta_x *= -np.sign(bond[2])
    Rx = rot_x(theta_x)
    xyz = np.array([np.dot(Rx,xyz_i) for xyz_i in xyz])
    
    return xyz

def read_geoms(direc,file,fix_ats=None):
    '''read in an xyz file with multiple geometries. Returns list of GEOM class instances'''
    f = open(direc+file)
    lines = f.readlines()
    f.close()
    
    Nat = int(lines[0])
    IMs = [] # list of image class instances
    
    ## read in data
    i=0
    while lines[i].strip().startswith(str(Nat)):
        try:
            E = float(lines[i+1].split()[-1])
        except:
            E = 0. # in case there is no energy
        atoms = []
        xyz = []
        for j in range(Nat):
            at, x, y, z = lines[i+2+j].split()
            atoms.append(at)
            xyz.append([float(x), float(y), float(z)])
        xyz = np.array(xyz)
        # orient to common system
        if fix_ats != None:
            xyz = orient(xyz,fix_ats)
        
        IMs.append(GEOM(atoms,xyz,E))
        i += Nat+2
        if i >= len(lines):
            break
        
    return IMs # returns list of GEOM class instances

#%% classes

class GEOM:
    '''Class for single geometry'''
    def __init__(self,atoms,xyz,E):
        '''atom labels (list of strings), coordinates (Nx3 dim np.array), energy (float)'''
        self.atoms = atoms
        self.masses = np.array([a2m[X] for X in self.atoms])
        self.xyz = xyz # Nx3 array
        self.E = E
        self.N_at = len(atoms)
        
        # if bonds are given, use them, otherwise create them via a proximity check
        bonds = []
        for a1 in range(self.N_at):
            for a2 in range(a1+1,self.N_at):
                p1 = self.xyz[a1]
                p2 = self.xyz[a2]
                if LA.norm(p2-p1) < 1.55: # covalent bonds
                    bonds.append([a1,a2])
        self.bonds = np.array(bonds)
    
    def as_3N(self):
        '''returns Nx3 xyz array as 3N array'''
        return np.hstack(self.xyz)
    
    def CoM(self):
        '''returns the center of mass of the conformer'''
        return np.sum(self.masses*self.xyz.T,axis=1)/np.sum(self.masses)
    
    def R(self,i,j):
        '''calculates distance between atoms i and j (counting from 0)'''
        return np.linalg.norm(self.xyz[i]-self.xyz[j])

    def A(self,i,j,k):
        '''calculates the angle formed by atoms i-j-k (counting from 0)'''
        v1 = self.xyz[i] - self.xyz[j]
        v2 = self.xyz[k] - self.xyz[j]
        return np.degrees(np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
    
    def D(self,i,j,k,l):
        '''calculates the dighedral formed by atoms i-j-k-l (counting from 0)'''
        # http://azevedolab.net/resources/dihedral_angle.pdf
        # three plane defining vectors
        q1 = self.xyz[i] - self.xyz[j]
        q2 = self.xyz[k] - self.xyz[j]
        q3 = self.xyz[l] - self.xyz[k]
        # normal vectors of two planes
        n1 = np.cross(q1,q2) / np.linalg.norm(np.cross(q1,q2))
        n2 = np.cross(-q2,q3) / np.linalg.norm(np.cross(-q2,q3))
        # orthogonal unity vectors
        u1 = n2
        u3 = q2/np.linalg.norm(q2)
        u2 = np.cross(u3,u1)
        # dihedral via theta=-atan2(n1*u2/n1*u1)
        return -np.degrees(np.math.atan2(np.dot(n1,u2),np.dot(n1,u1)))
    
    def dist_mat(self):
        '''returns distance matrix of xyz, i.e., a matrix where D_ij is the distance of atoms i and j'''
        D = np.array([[ LA.norm(xyz_i-xyz_j) for xyz_i in self.xyz] for xyz_j in self.xyz])
        return D
    
    def plot_geom(self,lab1=''):
        '''plots the xyz structure'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        
        # atoms
        for a in range(self.N_at):
            ax.scatter(self.xyz[a,0],self.xyz[a,1],self.xyz[a,2],color=at_cols[self.atoms[a]],
                       edgecolor='k',alpha=1.0,s=self.masses[a]*5)
        ax.set_box_aspect((np.ptp(self.xyz[:,0]), np.ptp(self.xyz[:,1]), np.ptp(self.xyz[:,2])))
        
        # bonds
        for a1,a2 in self.bonds:
            p1 = self.xyz[a1]
            p2 = self.xyz[a2]
            pm = p1+0.5*(p2-p1) # center between bond for color change
            ax.plot([p1[0],pm[0]],[p1[1],pm[1]],[p1[2],pm[2]],at_cols[self.atoms[a1]],alpha=0.8)
            ax.plot([pm[0],p2[0]],[pm[1],p2[1]],[pm[2],p2[2]],at_cols[self.atoms[a2]],alpha=0.8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title(lab1)
        plt.show()
    
    def write_ORCA_inp(self,direc,fname,method):
        '''uses the xyz file to create an ORCA input file with the given method'''
        f = open(direc+fname+'.inp','w+')
        for line in method:
            f.write(line + '\n')
        for a in range(self.N_at):
            f.write(' %s         %+.10f    %+.10f    %+.10f' %(self.atoms[a],*self.xyz[a]) + '\n')
        f.write('*' + '\n')
        f.close()


class ENS:
    '''class of the ensemble of geometries'''
    def __init__(self,direc,file,fix_ats=None):
        self.direc = direc
        self.file = file
        self.fix_ats = fix_ats # list of three atoms to reorient
        self.GEOMS = read_geoms(direc,file,fix_ats) # read in the geometries and (potentially) reorient
        self.Nconf = len(self.GEOMS) # number of conformers
        self.atoms = self.GEOMS[0].atoms # list of atom labels
        self.masses = self.GEOMS[0].masses # list of atom masses
        self.N_at = len(self.atoms) # number of atoms
        self.ener = np.array([X.E for X in self.GEOMS]) # energies in Eh
        self.E_rel = self.ener - np.min(self.ener) # relative energies in Eh
        
        # create directories
        if not os.path.isdir(direc+subfol):
            os.mkdir(self.direc+subfol)
    
    def write_xyz(self,C,lab):
        '''writes an ensemble xyz file with all conformers in cluster C (list of indices).'''
        file = open(self.direc+subfol+'/'+self.file[:-4] + '_' + lab + '.xyz', 'w+')
        for N in C:
            file.write('  %i' %self.N_at + '\n')
            file.write('        %+.8f' %self.ener[N] + '\n')
            xyz_N = self.GEOMS[N].xyz
            for a in range(self.N_at):
                file.write(' %s         %+.10f    %+.10f    %+.10f' %(self.atoms[a],*xyz_N[a]) + '\n')
        file.close()
    
    def write_ORCA_inp(self,direc,fpre,i_list,method):
        '''calls GEOM.write_ORCA_inp function for all conformers given in i_list'''
        fmt = '%0'+'%i'%(math.ceil(np.log10(self.Nconf))) + 'd'
        for i in i_list:
            self.GEOMS[i].write_ORCA_inp(direc+subfol+'/', fpre+'_C'+fmt %i+'_OptFreq', method)
    
    def distance_dist(self,i,j,Nbins=40):
        '''distribution of the angle between atoms i,j,k over the ensemble'''
        RR = np.array([X.R(i,j) for X in self.GEOMS])
        plt.figure()
        plt.hist(RR,bins=Nbins,density=True)
        plt.xlabel(r'distance between %s(%i)-%s(%i) / $\AA$' %(self.atoms[i],i,self.atoms[j],j))
        plt.ylabel('frequency')
        plt.yticks([])
        plt.show()

    def angle_dist(self,i,j,k,Nbins=40):
        '''distribution of the angle between atoms i,j,k over the ensemble'''
        AA = np.array([X.A(i,j,k) for X in self.GEOMS])
        plt.figure()
        plt.hist(AA,bins=Nbins,density=True)
        plt.xlabel('angle between %s(%i)-%s(%i)-%s(%i) / degree' %(self.atoms[i],i,self.atoms[j],j,self.atoms[k],k))
        plt.ylabel('frequency')
        plt.yticks([])
        plt.show()
    
    def dihed_dist(self,i,j,k,l):
        '''distribution of the angle between atoms i,j,k over the ensemble'''
        #https://stackoverflow.com/questions/22562364/circular-polar-histogram-in-python
        Ddeg = np.array([X.D(i,j,k,l) for X in self.GEOMS])
        Drad = np.deg2rad(Ddeg)

        Nbins = 80
        theta = np.linspace(-np.pi, np.pi, Nbins+1, endpoint=True)
        width = (2*np.pi) / Nbins
        base = 0.5
        radii = np.histogram(Drad,bins=theta,density=True)[0] # frequency
        radii *= (1.-base)/radii.max()
        
        plt.figure(figsize=(6,6))
        ax = plt.subplot(111, polar=True)
        ax.set_title('dihedral between %s(%i)-%s(%i)-%s(%i)-%s(%i)'
                  %(self.atoms[i],i,self.atoms[j],j,self.atoms[k],k,self.atoms[l],l))
        
        ax.bar(theta[:-1], radii, width=width, bottom=base, color='C0')
        ax.plot(Drad,np.ones(self.Nconf)*(base*0.8),'C0o', alpha=0.4)
        
        ax.set_yticks([0.,base,1.])
        ax.set_yticklabels(['','',''])
        ax.set_xticks(np.deg2rad([0.,60.,120.,180.,240.,300.]))
        ax.set_xticklabels(['0°','60°','120°',r'$\pm$180°',r'$-$120°',r'$-$60°'])
        plt.show()
    
    def ener_dist(self,Nbins=40):
        '''histogram of the energies'''
        AA = np.array([X.E for X in self.GEOMS])
        plt.figure()
        plt.hist(AA,bins=Nbins,density=True)
        plt.xlabel(r'energies / $E_h$')
        plt.ylabel('frequency')
        plt.yticks([])
        plt.show()
    
    def plot_geoms(self,i,j):
        '''plots two structures i and j'''
        CS_score = self.cos_sim_pair(i,j)
        DM_score = self.dist_mat_pair(i,j)
        xyz1 = self.GEOMS[i].xyz
        xyz2 = self.GEOMS[j].xyz
        lab1 = 'C%i, E=%.3f eV' %(i,self.E_rel[i]*Eh2eV)
        lab2 = 'C%i, E=%.3f eV' %(j,self.E_rel[j]*Eh2eV)
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title('cos similarity: %.1f' %(CS_score*100) + '%\n' +
                     r'dist mat distance: %.1f $\AA$' %(DM_score))
        
        # atoms
        for a in range(1,self.N_at):
            ax.scatter(xyz1[a,0],xyz1[a,1],xyz1[a,2],color='r',edgecolor='k',alpha=0.8,s=self.masses[a]*5)
            ax.scatter(xyz2[a,0],xyz2[a,1],xyz2[a,2],color='b',edgecolor='k',alpha=0.8,s=self.masses[a]*5)
        ax.set_box_aspect((np.ptp(xyz1[:,0]), np.ptp(xyz1[:,1]), np.ptp(xyz1[:,2])))
        ax.scatter(xyz1[0,0],xyz1[0,1],xyz1[0,2],color='r',edgecolor='k',alpha=0.8,
                   s=a2m[self.atoms[a]]*5, label=lab1)
        ax.scatter(xyz2[0,0],xyz2[0,1],xyz2[0,2],color='b',edgecolor='k',alpha=0.8,
                   s=a2m[self.atoms[a]]*5, label=lab2)
        
        # bonds
        for a1,a2 in self.GEOMS[i].bonds:
            p11 = xyz1[a1]
            p21 = xyz1[a2]
            ax.plot([p11[0],p21[0]],[p11[1],p21[1]],[p11[2],p21[2]],'r',alpha=0.8)
        for a1,a2 in self.GEOMS[j].bonds:
            p12 = xyz2[a1]
            p22 = xyz2[a2]
            ax.plot([p12[0],p22[0]],[p12[1],p22[1]],[p12[2],p22[2]],'b',alpha=0.8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        legend_elements=[Line2D([0], [0], marker='o', color='w', label=lab1,
                              markerfacecolor='r', markersize=10),
                         Line2D([0], [0], marker='o', color='w', label=lab2,
                              markerfacecolor='b', markersize=10)]
        ax.legend(handles=legend_elements, loc='upper right')
        limits = np.array([getattr(ax, f'get_{axis}lim')() for axis in 'xyz'])
        ax.set_box_aspect(np.ptp(limits, axis = 1))
        plt.show()

    def plot_geoms_many(self,ilist):
        '''plots of all geoms in ilist'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
 
        Xi = [self.GEOMS[i] for i in ilist]
        alph = 1./np.sqrt(len(Xi))
        for X in Xi:
            xyz1 = X.xyz
            
            # atoms
            for a in range(self.N_at):
                ax.scatter(xyz1[a,0],xyz1[a,1],xyz1[a,2],color=at_cols[self.atoms[a]],
                           edgecolor='k',alpha=alph,s=self.masses[a]*5)
            # bonds
            for a1,a2 in X.bonds:
                p1 = xyz1[a1]
                p2 = xyz1[a2]
                pm = p1+0.5*(p2-p1) # half way between atoms a1 and a2
                ax.plot([p1[0],pm[0]],[p1[1],pm[1]],[p1[2],pm[2]],at_cols[self.atoms[a1]],alpha=alph*0.8)
                ax.plot([pm[0],p2[0]],[pm[1],p2[1]],[pm[2],p2[2]],at_cols[self.atoms[a2]],alpha=alph*0.8)
            
        #ax.set_box_aspect((np.ptp(xyz1[:,0]), np.ptp(xyz1[:,1]), np.ptp(xyz1[:,2])))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.view_init(elev=0., azim=0.) # look along aligned atoms
        limits = np.array([getattr(ax, f'get_{axis}lim')() for axis in 'xyz'])
        ax.set_box_aspect(np.ptp(limits, axis = 1))
        plt.show()
    
    def cos_sim_pair(self,i,j):
        '''calculates the cosine similarity for conformers i and j'''
        R1 = self.GEOMS[i].masses*LA.norm(self.GEOMS[i].xyz-self.GEOMS[i].CoM(),axis=1)
        R2 = self.GEOMS[j].masses*LA.norm(self.GEOMS[j].xyz-self.GEOMS[j].CoM(),axis=1)
        R1 = np.sort(R1)
        R2 = np.sort(R2)
        s = np.dot(R1,R2)/(LA.norm(R1)*LA.norm(R2))
        if s > 1.:
            s = 1.
        elif s < -1.:
            s = -1.
        return 1.-np.arccos(s)/np.pi
    
    def cos_sim_all(self,plotting=True):
        '''calculates the cosine similarity matrix for all the structures'''
        #CoMs = np.array([ X.CoM() for X in self.GEOMS ])
        Rs = np.array([X.masses*LA.norm(X.xyz-X.CoM(),axis=1) for X in self.GEOMS])
        normRs = LA.norm(Rs, axis=1)
        
        #s_mat = np.array([[np.dot(R1,R2)/(LA.norm(R1)*LA.norm(R2)) for R1 in Rs] for R2 in Rs])
        s_mat = np.identity(self.Nconf)
        for i in range(self.Nconf):
            for j in range(i+1,self.Nconf):
                s = np.dot(Rs[i],Rs[j])/(normRs[i]*normRs[j])
                if s > 1.:
                    s = 1.
                elif s < -1.:
                    s = -1.
                s_mat[i,j] = s
                s_mat[j,i] = s
        score = 1.-np.arccos(s_mat)/np.pi
        
        if plotting:
            fig, ax = plt.subplots(figsize=(7,6))
            cax = ax.matshow(score*100)
            fig.colorbar(cax,label='cosine similarity (%)')
            plt.show()
        
        return score
    
    def dist_mat_pair(self,i,j):
        '''calculates the similarity of the distance matrices (GEOM.dist_mat()) of conformers i and j'''
        D1 = self.GEOMS[i].dist_mat() # diagonal is zero
        D2 = self.GEOMS[j].dist_mat() # diagonal is zero
        dD = np.abs(D1-D2) # diagonal is zero
        score = np.sum( dD[np.triu_indices(self.N_at)] ) # similarity score
        return score
    
    def dist_mat_all(self,plotting=True):
        '''calculates the similarity of the distance matrices (GEOM.dist_mat()) for all conformers'''
        Ds = [X.dist_mat() for X in self.GEOMS]
        score = np.zeros((self.Nconf,self.Nconf))
        for i in range(self.Nconf):
            for j in range(i+1,self.Nconf):
                s = np.sum( np.abs(Ds[i]-Ds[j])[np.triu_indices(self.N_at)] )
                score[i,j] = s
        score = score + score.T # score[i,j] = score[j,i]
        
        if plotting:
            fig, ax = plt.subplots(figsize=(7,6))
            cax = ax.matshow(score)
            fig.colorbar(cax,label='distance similarity score')
            plt.show()
        
        return score
    
    def plot_dendrogram(self,which='dist_mat',method='single'):
        '''plotting dendrogram using scipy.cluster.hierarchy using either self.dist_mat_all (which='dist_mat'),
        self.cos_sim_all (which='cos_sim') similarity scores or just the energy GEOM.E (which='energy').
        <method> should be one of the ones reported in scipy.cluster.hierarchy.linkage.'''
        if which == 'dist_mat':
            X = self.dist_mat_all(plotting=False)
            ylab = 'distance matrix similarity'
        elif which == 'cos_sim':
            X = (1.-self.cos_sim_all(plotting=False))*100.
            ylab = 'cosine similarity distance (%)'
        elif which == 'energy':
            X = np.array([np.abs(self.E_rel-self.E_rel[i]) for i in range(self.Nconf)])
            ylab = 'rel. energy'
        else:
            raise ValueError('<which> needs to be "dist_mat", "cos_sim", or "energy"!')
            
        Xcondensed = X[np.triu_indices(self.Nconf,1)] # only upper triangle w/o diagonal
        res = hierarchy.linkage(Xcondensed,method=method)
        plt.figure(figsize=(15,5))
        hierarchy.dendrogram(res)
        plt.xlabel('index of conformer')
        plt.ylabel(ylab)
        ydat = res[:,2]
        yrange = np.max(ydat) - np.min(ydat)
        plt.ylim(np.max([0,np.min(ydat)-yrange*0.05]), np.max(ydat)+yrange*0.05)
        plt.tight_layout()
        plt.savefig(self.direc+subfol+'/'+'dendrogram_%s_%s.png' %(which,method),
                    format='png',dpi=300,bbox_inches='tight')
        plt.show()
        
        return res
    
    def plot_geoms_cluster(self,dn,N,t):
        '''uses self.plot_geoms_many() to plot all structures in the same cluster as conformer N.
        t is the threshold value what to consider a cluster (y-value in the dendrogram).
        Larger t will increase the number of conformers in the cluster of conformer N'''
        if N > self.Nconf:
            raise ValueError('N needs to be one of the %i conformer indices!' %self.Nconf)
        
        T = hierarchy.fcluster(dn, t, criterion='distance') # flatten clusters
        cluster = T[N] # which cluster is conformer N in
        which = np.arange(self.Nconf)[T==cluster] # find all conformers that are in cluster
        print('Plotting cluster %i/%i with %i geometries.' %(cluster, np.max(T), len(which)))
        print('I.e., conformers ',  which)
        self.plot_geoms_many(which)
        return T
    
    def reduce_ENS(self,dn,t):
        '''reduces the structures of the ensemble according to the clustering. All structures
        within a cluster (as defined by t) will be represented by their lowest energy member.'''
        T = hierarchy.fcluster(dn, t, criterion='distance') # flatten clusters
        Ncl = np.max(T) # will correspond to the resulting number of structures
        repres = [] # list of representatives (lenght of Ncl)
        for iCl in range(1,Ncl+1):
            which = np.arange(self.Nconf)[T==iCl] # find all conformers that are in cluster
            rep_i = which[ np.argmin( [self.GEOMS[i].E for i in which] ) ] # lowest energy representative
            print('representing cluster %02d through conformer %03d' %(iCl,rep_i))
            repres.append(rep_i)
        self.write_ORCA_inp(direc,'1-NG_rep',repres,method) # method is defined at the top
        print('\nORCA inp files written!')


######################################################################
direc = os.getcwd()+'/'
ENS_conf = ENS(direc,'crest_conformers.xyz',fix_ats=[2,11,0])
#ENS_rot  = ENS(direc,'crest_rotamers.xyz',fix_ats=[2,11,0])

#%% write reordered structures (much easier to visualize with Avogadro)
ENS_conf.write_xyz(np.arange(ENS_conf.Nconf),'reorient')
#ENS_rot.write_xyz( np.arange(ENS_rot.Nconf), 'reorient')

#%% show distributions of energies or internal coordinates
#ENS_conf.ener_dist()
#ENS_conf.distance_dist(15,10)
#ENS_conf.angle_dist(9,8,0)
ENS_conf.dihed_dist(14,5,3,12)

#%% 3D plotting of geometries
ENS_conf.plot_geoms(0,1) # overlay two structures, colored
ENS_conf.plot_geoms_many([0,1,2,3]) # list of structures (list can also be just one element)

#%% write ORCA inp for a list of conformers
#ENS_conf.write_ORCA_inp(direc,'1-NG',[0,2,3,4,5,7,10],method) # method is defined at the top

#%% comparison metrics
#ENS_conf.cos_sim_pair(0,1)
ENS_conf.cos_sim_all()
#ENS_conf.dist_mat_pair(0,1)
#ENS_conf.dist_mat_all()

#%% hierarchy clustering ##
dn = ENS_conf.plot_dendrogram(which='cos_sim',method='weighted')

#%% automatic reduction according to the dendogram clustering
t = 7.4 # distance threshold under which to cluster (y-axis of dendogram); for cos_sim, t=2 means 98% cos_sim
ENS_conf.reduce_ENS(dn,t)

#%% investigate similarity manually within one specific cluster (read t and iat from the dendogram)
t = 7.4 # distance threshold under which to cluster (y-axis of dendogram)
iat = 0 # any conformer inside the cluster
# 3D plot of all structures in that cluster (overlayed)
T = ENS_conf.plot_geoms_cluster(dn, iat, t)
# or alternatively: write one xyz file for all structures in that cluster for investigation via Avogadro
iC = T[iat] # cluster number
ENS_conf.write_xyz(np.arange(ENS_conf.Nconf)[T==iC], 't%i_C%i' %(t,iC))
# --> this could be used for manual reduction (take one strucutre from each cluster as representative)