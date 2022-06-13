import matplotlib.pyplot as plt
import numpy as np
import dcd as dcd
import sys
 

if len(sys.argv)<2:
	print("Usage --> python3 3d_quarter_distribution.py <input dcd> <index of center of mass of ligand> <no of replicas> <dwell_x> <dwell_y> <dwell_Z> <r1> <r2>")
if len(sys.argv)==9:
	input_dcd= sys.argv[1]
	com_index=int(sys.argv[2])
	n_replicas = int(sys.argv[3])
	dwell_x=float(sys.argv[4])
	dwell_y=float(sys.argv[5])
	dwell_z=float(sys.argv[6])
	r1=float(sys.argv[7])
	r2=float(sys.argv[8])

	cc,co=dcd.read_dcd_data_ref(input_dcd)
	print ("Center of mass index	="+str(com_index))
	print ("no. of replicas		="+str(n_replicas))
	n_atm_per_res=int(int(np.size(co,1))/n_replicas)
	
	print("no. of atoms per residue	="+str(n_atm_per_res))
	print("no. of frames		="+str(np.size(co,0)))
	print("      ")
	print(np.shape(co))
	n_frames=int(int(np.size(co,0)))
	co=co[:, com_index:-1:n_atm_per_res, :]
	print(np.shape(co))
	#saving every 5th frame from first 1000 frames

#	co=co[50:2000:2, :, :]      #first 500ns frames
	co=co[-2000:-1:2, :, :]     #last 500ns frames
	print (np.shape(co))
	points=0
	q1=0
	q2=0
	q3=0
	q4=0
	for p in range(0,200):
		for q in range(0,100):
			x=float(co[p][q][0])
			y=float(co[p][q][1])
			z=float(co[p][q][2])
			dx=x-dwell_x
			dy=y-dwell_y
			dz=z-dwell_z
			
			r=((dx*dx)+(dy*dy)+(dz*dz))**0.5
			if r>r1 and r<r2:
				points=points+1
				if dx>0 and dy>0:
					q1=q1+1				
				elif dx>0 and dy<0:
					q2=q2+1
				elif dx<0 and dy<0:
					q3=q3+1
				elif dx<0 and dy>0:
					q4=q4+1

q1=(q1/points)*100
q2=(q2/points)*100
q3=(q3/points)*100
q4=(q4/points)*100

print("points"+str(points))
print ("q1 --> "+ str(q1))
print ("q2 --> "+ str(q2))
print ("q3 --> "+ str(q3))
print ("q4 --> "+ str(q4))
