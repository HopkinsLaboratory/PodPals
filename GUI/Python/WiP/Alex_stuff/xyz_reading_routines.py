
a2m = {'H' :  1.,                                                              'He':  4.,
       'Li':  7., 'Be':  9., 'B': 11, 'C': 12., 'N': 14., 'O': 16., 'F' : 19., 'Ne': 20.,
       'Na': 23., 'Mg': 24.,                    'P': 31., 'S': 32., 'Cl': 35., 'Ar': 40.,
       } # from atom label to mass

at_cols = {'C': 'black', 'H': 'dimgrey', 'O': 'orangered', 'N': 'royalblue',
           'S': 'gold', 'F': 'cyan', 'Cl': 'lime', 'Br': 'firebrick'}

def parse(file,off=5):
    '''Parse the molecular geometries from a file (.inp, .gjf, .xyz).
    Set off to the number of lines before the geometry starts (e.g. off=5 for .gjf)'''
    f = open(file)
    lines = f.readlines()
    f.close()
    
    labs = []
    mass = []
    xyz = []
    for i in range(off,len(lines)):
        if len(lines[i].split()) == 4:
            lab, x, y, z = lines[i].split()
            mass.append(a2m[lab])
            labs.append(lab)
            xyz.append([float(x), float(y), float(z)])
    
    return labs, np.array(mass), np.array(xyz)

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
    '''Orientes the molecule such that the fix_ats (array of length 3) have the following
    properties: 1st one lies at (0,0,0), 2nd one lies at (+x,0,0), 3rd one lies at (x,+y,0)'''
    xyz = np.copy(xyz_in)
    ## translate fix_ats[0] to origin
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


def plot_1geom(xyz1,lab1,cov_thresh=1.55):
    '''plots two structures xyz1 and xyz2'''
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    
    # atoms
    for a in range(N_at):
        ax.scatter(atoms[0,a],atoms[1,a],atoms[2,a],color=at_cols[atom_labs[a]],
                   edgecolor='k',alpha=1.0,s=mass[a]*5)
    ax.set_box_aspect((np.ptp(xyz1[:,0]), np.ptp(xyz1[:,1]), np.ptp(xyz1[:,2])))
    
    # bonds
    #bonds = parse_bonds()
    for a1 in range(N_at):
        for a2 in range(a1+1,N_at):
            p1 = xyz1[a1]
            p2 = xyz1[a2]
            pm = p1+0.5*(p2-p1) # half way between atoms a1 and a2
            if LA.norm(p21-p11) < cov_thresh: # covalent bonds
                ax.plot([p1[0],pm[0]],[p1[1],pm[1]],[p1[2],pm[2]],at_cols[atom_labs[a1]],alpha=0.8)
            	ax.plot([pm[0],p2[0]],[pm[1],p2[1]],[pm[2],p2[2]],at_cols[atom_labs[a2]],alpha=0.8)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.title(lab1)
    plt.show()

def plot_2geoms(xyz1,xyz2,lab1,lab2,cov_thresh=1.55):
    '''plots two structures xyz1 and xyz2'''
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    
    # atoms
    for a in range(1,N_at):
        ax.scatter(xyz1[a,0],xyz1[a,1],xyz1[a,2],color='r',edgecolor='k',alpha=0.8,s=mass[a]*5)
        ax.scatter(xyz2[a,0],xyz2[a,1],xyz2[a,2],color='b',edgecolor='k',alpha=0.8,s=mass[a]*5)
    ax.set_box_aspect((np.ptp(xyz1[:,0]), np.ptp(xyz1[:,1]), np.ptp(xyz1[:,2])))
    ax.scatter(xyz1[0,0],xyz1[0,1],xyz1[0,2],color='r',edgecolor='k',alpha=0.8,
               s=mass[0]*5, label=lab1)
    ax.scatter(xyz2[0,0],xyz2[0,1],xyz2[0,2],color='b',edgecolor='k',alpha=0.8,
               s=mass[0]*5, label=lab2)
    
    # bonds
    #bonds = parse_bonds()
    for a1 in range(N_at):
        for a2 in range(a1+1,N_at):
            p11 = xyz1[a1]
            p21 = xyz1[a2]
            p12 = xyz2[a1]
            p22 = xyz2[a2]
            if LA.norm(p21-p11) < cov_thresh: # covalent bonds
                ax.plot([p11[0],p21[0]],[p11[1],p21[1]],[p11[2],p21[2]],'r',alpha=0.8)
            if LA.norm(p22-p12) < cov_thresh: # covalent bonds
                ax.plot([p12[0],p22[0]],[p12[1],p22[1]],[p12[2],p22[2]],'b',alpha=0.8)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    legend_elements=[Line2D([0], [0], marker='o', color='w', label=lab1,
                          markerfacecolor='r', markersize=10),
                     Line2D([0], [0], marker='o', color='w', label=lab2,
                          markerfacecolor='b', markersize=10)]
    ax.legend(handles=legend_elements, loc='upper right')
    plt.show()
