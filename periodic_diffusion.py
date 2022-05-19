# save pdb of 152 frames
#pdb file of trajectory must start with lines having coordinates, i.e., no extra headers
#delete pdb's header before
#for 1 atom pdb only

import sys
import numpy as np
from vectors import Point, Vector
import itertools

MSD = np.array([])
SUM_dr_2 = 0.
if len(sys.argv)< 2:
    print ("usage: python periodic_diffusion.py <input_file> <no. of frames>")

if len(sys.argv)==3:	
    pdb = sys.argv[1]
    n_frames =sys.argv[2]
    n_frames =(int(n_frames)-1) 
#    print (n_frames)
    n_ligands = 1
    n_atm_ligand = 1
#    print ("input = "+ pdb + ", " "no. of ligands = "+ str(n_ligands) + " and "+ "no. of atoms in each ligand = "+ str(n_atm_ligand))
    END = (n_ligands * n_atm_ligand)+1.
    f1=open(pdb)

    for i, l in enumerate (f1):
    	if l.startswith('ATOM'):
            x_1 = float(l[30:38])
            y_1 = float(l[38:46])
            z_1 = float(l[46:54])
            p_one = Point(x_1, y_1, z_1)

            f2 = open(pdb)
            for i2, l2 in enumerate (f2):
                if i2 == i+END:					
                    x_2 = float(l2[30:38])
                    y_2 = float(l2[38:46])
                    z_2 = float(l2[46:54])
                    p_two = Point(x_2, y_2, z_2)
                    dx = x_1 - x_2
                    dy = y_1 - y_2
                    dz = z_1 - z_2
#edit 400. according to size of periodic box					
                    if dx < -300.:
                        dx = dx + 400.  
                    if dy < -300.:
                        dy = dy + 400.   
                    if dz < -300.:
                        dz = dz + 400. 
                    if dx > 300.:
                        dx = dx - 400.	
                    if dy > 300.:
                        dy = dy -400.	
                    if dz > 300.:
                        dz = dz - 400.
			
                    dr_2 = dx*dx + dy*dy + dz*dz
                    MSD = np.append(MSD, dr_2)	


    MSD = MSD.reshape(-1, 1)      
    MSD = np.mean(MSD, axis=1)
    MSD = MSD.reshape(-1, 1)        
    MSD = np.cumsum(MSD, axis = 0)
#    print (MSD)
    SUM_MSD = MSD[-1]
    msd = np.mean(SUM_MSD)

    D = msd/((2.*3.)*(0.05*10000.*n_frames))
    print (str(D/10000))






