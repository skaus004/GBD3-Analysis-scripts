import MDAnalysis as mda
from MDAnalysis import *

import MDAnalysis.coordinates
import MDAnalysis.topology.TOPParser
import numpy.linalg
import sys
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import distances

#constant_values
kB=1.38e-23  #m^2Kg/s^2.K
T=298 #K
n=0.0008921 #Kg/m.s

ligand = sys.argv[1]

u = mda.Universe(ligand)  #open pdb
print ("Rgyr: {0:g} A".format(u.atoms.radius_of_gyration()))
r= u.atoms.radius_of_gyration()
centroid = u.atoms.center_of_geometry()

Diff_Coeff=(kB*T)/(6e-14*3.14*n*r)
print("Diffusion Coefficient ="+str(Diff_Coeff)+" cm^2/s")
print("Center of molecule =" + str(centroid))
