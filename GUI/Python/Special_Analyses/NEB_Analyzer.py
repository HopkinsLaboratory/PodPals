# -*- coding: utf-8 -*-
"""
Script to read in ORCA NEB calculation and plot data.

Created on Wed Jan 26 17:07:00 2022

@author: Alexander Haack
"""

import os
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

a2m = {'H' :  1.,                                                              'He':  4.,
       'Li':  7., 'Be':  9., 'B': 11, 'C': 12., 'N': 14., 'O': 16., 'F' : 19., 'Ne': 20.,
       'Na': 23., 'Mg': 24.,                    'P': 31., 'S': 32., 'Cl': 35., 'Ar': 40.,
       } # from atom label to mass

at_cols = {'C': 'black', 'H': 'dimgrey', 'O': 'orangered', 'N': 'royalblue',
           'S': 'gold', 'F': 'cyan', 'Cl': 'lime', 'Br': 'firebrick', 'P': 'orange'}

def dist(xyz1,xyz2):
    '''calcualtes the distance between two 3N vectors'''
    return LA.norm(xyz1-xyz2)

def read_path(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    Nat = int(lines[0])
    IMs = [] # list of image class instances
    
    ## read in data
    i=0
    while lines[i].startswith(str(Nat)):
        E = float(lines[i+1].split()[-1])
        atoms = []
        xyz = []
        for j in range(Nat):
            at, x, y, z = lines[i+2+j].split()
            atoms.append(at)
            xyz.append([float(x), float(y), float(z)])
        IMs.append(Imag(atoms,np.array(xyz),E))
        i += Nat+2
        if i >= len(lines):
            break
        
    return path(IMs) # return path class instance

class Imag:
    '''Class for single image of the path'''
    def __init__(self,atoms,xyz,E,bond_thresh=1.518):
        '''atom labels (list of strings), coordinates (Nx3 dim np.array), energy (float)'''
        self.atoms = atoms
        self.xyz = xyz
        self.E = E
        self.N_at = len(atoms)
        
        bonds = []
        for a1 in range(self.N_at):
            for a2 in range(a1+1,self.N_at):
                p1 = self.xyz[a1]
                p2 = self.xyz[a2]
                if LA.norm(p2-p1) < bond_thresh: # covalent bonds
                    bonds.append([a1,a2])
        self.bonds = np.array(bonds)
    
    def as_3N(self):
        '''returns Nx3 xyz array as 3N array'''
        return np.hstack(self.xyz)
    
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
    
    def export_xyz(self,name):
        '''exports the coordinates of this image to an xyz file'''
        f = open(name,'w+')
        f.write('%i' %(self.N_at) + '\n')
        f.write('Coordinates of Image, E %.8f ' %(self.E)  + '\n')
        for i in range(self.N_at):
            f.write('  %s  %+9.6f,  %+9.6f,  %+9.6f' %(self.atoms[i],*(self.xyz[i]),) + '\n')
        f.close()
    
    def plot_geom(self,lab=''):
        '''plots the xyz structure'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        
        # atoms
        for a in range(self.N_at):
            try:
                col = at_cols[self.atoms[a]] # coloring of the atom
            except:
                col = 'red'
            try:
                size = a2m[self.atoms[a]]*5
            except:
                size = 60
            ax.scatter(self.xyz[a,0],self.xyz[a,1],self.xyz[a,2],color=col,
                       edgecolor='k',alpha=1.0,s=size)
        ax.set_box_aspect((np.ptp(self.xyz[:,0]), np.ptp(self.xyz[:,1]), np.ptp(self.xyz[:,2])))
        
        # bonds
        for a1,a2 in self.bonds:
            p1 = self.xyz[a1]
            p2 = self.xyz[a2]
            pm = p1+0.5*(p2-p1) # center between bond for color change
            try:
                col1 = at_cols[self.atoms[a1]]
                col2 = at_cols[self.atoms[a2]]
            except:
                col1 = 'k'
                col2 = 'k'
            ax.plot([p1[0],pm[0]],[p1[1],pm[1]],[p1[2],pm[2]],col1,alpha=0.8)
            ax.plot([pm[0],p2[0]],[pm[1],p2[1]],[pm[2],p2[2]],col2,alpha=0.8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.title(lab)
        plt.show()

## end of image class

class path:
    '''Class of the MEP'''
    def __init__(self,IMs):
        self.IMs = IMs
        self.Nim = len(self.IMs)
        E = np.array([x.E for x in self.IMs])
        self.E = E-E.min()
        self.d = np.array([dist(self.IMs[0].as_3N(), self.IMs[i].as_3N()) for i in range(self.Nim) ])
        
        self.HEI = np.argmax(self.E)
        self.atoms = self.IMs[0].atoms
    
    def lab(self,*ind):
        llist = []
        for i in ind:
            llist.append(self.atoms[i])
            llist.append('$_{%i}$' %i)
            llist.append(',')
        label = '%s%s%s'*len(ind) %(*llist,)
        return label[:-1]
    
    def get_R(self,i,j):
        return np.array([IM.R(i,j) for IM in self.IMs])

    def get_A(self,i,j,k):
        return np.array([IM.A(i,j,k) for IM in self.IMs])

    def get_D(self,i,j,k,l,modulo=True):
        if modulo:
            D_path = np.array([IM.D(i,j,k,l)%360 for IM in self.IMs])
        else:
            D_path = np.array([IM.D(i,j,k,l) for IM in self.IMs])
        return D_path
    
    def export_xyz(self,i):
        '''export image i of the path as xyz file'''
        self.IMs[i].export_xyz(direc + file.replace('trj','Im%02d' %i))
    
    def plot_geom(self,i):
        '''creates 3D plot of the geometry of image i'''
        self.IMs[i].plot_geom('Image %i' %i)
    
    def plot(self,X,Y,xlab,ylab):
        '''plots X versus Y with respective labels'''
        if len(X) != len(Y):
            raise ValueError('X and Y arrays are not matching!')
        f = plt.figure(figsize=(6,6))
        f.tight_layout()
        plt.plot(X, Y, 'kx-', label='MEP')
        plt.plot(X[self.HEI], Y[self.HEI], 'rx', ms=15, label='HEI (Im_%i)' %self.HEI)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.legend()
        plt.show()
    
    def plot_EvDist(self):
        '''Plot energy versus distance between images'''
        self.plot(self.d,self.E,r'distance / $\AA$',r'rel. energy / $E_h$')

    def plot_EvR(self,i,j):
        '''Plot energy versus R(i,j) (counting from 0)'''
        self.plot(self.get_R(i,j),self.E,
                  r'R(%s) / $\AA$' %self.lab(i,j),r'rel. energy / $E_h$')
        
    def plot_EvA(self,i,j,k):
        '''Plot energy versus A(i,j,k) (counting from 0)'''
        self.plot(self.get_A(i,j,k),self.E,
                  r'A(%s) / degrees' %self.lab(i,j,k),r'rel. energy / $E_h$')
        
    def plot_EvD(self,i,j,k,l,modulo=True):
        '''Plot energy versus D(i,j,k,l) (counting from 0)'''
        self.plot(self.get_D(i,j,k,l,modulo=modulo),self.E,
                  r'D(%s) / degrees' %self.lab(i,j,k,l),r'rel. energy / $E_h$')
    
    def plot_RvDist(self,i,j):
        '''plot R(i,j) vs path progression (counting from 0)'''
        self.plot(self.d,self.get_R(i,j),r'distance / $\AA$',r'R(%s) / $\AA$' %self.lab(i,j))
    
    def plot_AvDist(self,i,j,k):
        '''plot A(i,j,k) vs path progression (counting from 0)'''
        self.plot(self.d,self.get_A(i,j,k),
                  r'distance / $\AA$',r'A(%s) / degrees' %self.lab(i,j,k))
    
    def plot_DvDist(self,i,j,k,l,modulo=True):
        '''plot D(i,j,k,l) vs path progression (counting from 0)'''
        self.plot(self.d,self.get_D(i,j,k,l,modulo=modulo),
                  r'distance / $\AA$',r'D(%s) / degrees' %self.lab(i,j,k,l))

## end of path class

class path_history:
    def __init__(self,direc,file,NIm):
        '''Specify the directory <direc> (str), the file name <file> (str),
        and the number of images in the NEB <NIm> (int)'''
        self.NIm = int(NIm)
        all_path = read_path(direc+file)
        NIm_tot = all_path.Nim
        if NIm_tot%NIm != 0:
            raise ValueError('Error, total number of images is not divisible by NIm!')
        else:
            self.Npaths = NIm_tot//self.NIm
            print('Found %i paths with %i images each!' %(self.Npaths,self.NIm))
        self.path_hist = [path([all_path.IMs[j] for j in range(i*NIm,(i+1)*NIm)]) for i in range(self.Npaths)]
    
    def plot_EvDist(self):
        ao = 0.2
        f = plt.figure(figsize=(6,6))
        f.tight_layout()
        for i,path_i in enumerate(self.path_hist):
            plt.plot(path_i.d, path_i.E, 'kx-', alpha=ao+(1-ao)*i/self.Npaths)
            plt.plot(path_i.d[path_i.HEI], path_i.E[path_i.HEI], 'rx', ms=15, alpha=ao+(1-ao)*i/self.Npaths)
        plt.xlabel(r'distance / $\AA$')
        plt.ylabel(r'rel. energy / $E_h$')
        #plt.legend()
        plt.show()
    
    def plot_EvDist_3D(self):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        fig.tight_layout()
        XX = np.array([ path_i.d for path_i in self.path_hist ])[::-1]
        YY = np.array([ np.arange(self.Npaths) for _ in range(self.NIm) ]).T
        ZZ = np.array([ path_i.E for path_i in self.path_hist ])[::-1]
        
        ax.plot_surface(XX,YY,ZZ, cmap='viridis')
        
        ax.set_xlabel(r'distance / $\AA$')
        ax.set_ylabel('MEP iteration')
        ax.set_zlabel(r'rel. energy / $E_h$')
        ax.set_yticks(np.arange(self.Npaths), ['%i'%i for i in range(self.Npaths,0,-1)])
        plt.show()
        

#######################################################################################
##  You only need to edit stuff here!                                                ##
##  Define directory and filename of the name_MEP_trj.xyz file                       ##
##  The name_MEP_trj.xyz file is created after each NEB iteration.                   ##
##  The calculation does not need to be finished to use this script!                 ##
##                                                                                   ##
##  At the bottom, you can pre-define plots:                                         ##
##   MEP.plot_EvDist()        - plots energy vs. path length                         ##
##   MEP.plot_EvR(i,j)        - plots energy vs. bond distance between atoms i-j     ##
##   MEP.plot_EvA(i,j,k)      - plots energy vs. angle between atoms i-j-k           ##
##   MEP.plot_EvD(i,j,k,l)    - plots energy vs. dihedral between atoms i-j-k-l      ##
##   MEP.plot_RvDist(i,j)     - plots bond distance between atoms i-j vs. path length##
##   MEP.plot_AvDist(i,j,k)   - plots angle between atoms i-j-k vs. path length      ##
##   MEP.plot_DvDist(i,j,k,l) - plots dihedral between atoms i-j-k-l vs. path length ##
##                                                                                   ##
##  All atom numbers i,j,k,l are counted starting from 0                             ##
##  You can also type these commands into the console.                               ##
##  Plots will show the path as well as the highest energy image (HEI)               ##
##                                                                                   ##
##  The MEP class also provides other functions. You can get the plotted arrays with ##
##   MEP.E              - returns an array with E along the path                     ##
##   MEP.get_R(i,j)     - returns an array with R(i,j) along the path                ##
##   MEP.get_A(i,j,k)   - returns an array with A(i,j,k) along the path              ##
##   MEP.get_D(i,j,k,l) - returns an array with D(i,j,k,l) along the path            ##
##                                                                                   ##
##  The images of the MEP can be exported as xyz via                                 ##
##   MEP.export_xyz(i)  - where i is the number of the image (starting from 0)       ##
##                      - in particular, i can be set to MEP.HEI for the TS guess    ##
##                                                                                   ##
##  The images of the MEP can be viewed as 3D plot via                               ##
##   MEP.plot_geom(i)  - where i is the number of the image (starting from 0)        ##
##                     - in particular, i can be set to MEP.HEI for the TS guess     ##
##                                                                                   ##
#######################################################################################

direc = os.getcwd()+'/'
file = 'ACE-H_mixed_1_1_a2d_NEB-TS_MEP_trj.xyz' # _MEP_trj.xyz file

MEP = read_path(direc+file) # don't change this

## plotting of path in different projections
MEP.plot_EvDist()
#MEP.plot_EvR(13,16)
#MEP.plot_EvA(14,15,21)
#MEP.plot_EvD(13,10,14,16,modulo=True)
#MEP.plot_DvDist(13,10,14,16,modulo=True)

#%% path history class for plotting the evolution of the NEB
file_all = 'ACE-H_mixed_1_1_a2d_NEB-TS_MEP_ALL_trj.xyz'  # _MEP_ALL_trj.xyz file
MEP_hist = path_history(direc,file_all,18) # initialize class instance
MEP_hist.plot_EvDist()
MEP_hist.plot_EvDist_3D()