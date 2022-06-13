import dcd as dcd
import numpy as np
import sys
import math

if len(sys.argv)<2:
        print("Usage --> python3 count_termination_dcd.py <input dcd> <index of center of mass of ligand> <no of replicas> <x-coordinate of dwelling sphere> <y...> <z...> <dwell sphere radius>")


if len(sys.argv)==8:
	input_dcd = sys.argv[1]
	com_index =int(sys.argv[2])
	n_replicas = int(sys.argv[3])
	d1 = int(float(sys.argv[4]))
	d2 = int(float(sys.argv[5]))
	d3 = int(float(sys.argv[6]))
	dwell_cut_off = int(float(sys.argv[7]))

	print ("Center of mass index     = "+str(com_index))
	print ("no. of replicas          = "+ str(n_replicas))
	cc,co=dcd.read_dcd_data_ref(input_dcd)  #co are the cordinates in numpy array [frames, atoms, xyz]==[axis0, axis1, axis2]
	n_atm_per_res=int(int(np.size(co, 1))/n_replicas)

	print ("no. of atoms per residue = "+ str(n_atm_per_res))
	print ("no. of frames            = "+str(np.size(co,0)))
	print ("     ")
	print (np.shape(co))
	n_frames=int(int(np.size(co,0)))
	co=co[:, com_index:-1:n_atm_per_res, :]  # [frames, com_atom, xyz] == [axis0, axis1, axis2]
	print (np.shape(co))

	####saving every 10th frame from last 10000 frames, so 1000 frames; skip rest
	co=co[-10000:-1:10, :, :]
	print (np.shape(co))
	
	n50=0
	n75=0
	n100=0
	n125=0
	n150=0
	n175=0
	n200=0
	n225=0
	n250=0
	n300=0
	for p in range(0,np.size(co,0)):
		for q in range(0,np.size(co,1)):
			x = int(co[p][q][0])
			y = int(co[p][q][1])
			z = int(co[p][q][2])
			r = (((x-d1)**2)+((y-d2)**2)+((z-d3)**2))**0.5
			if r < 50:
				n50=n50+1		
			if r < 75 and r>50:
				n75=n75+1 
			if r<100 and r>75:
				n100=n100+1
			if r<125 and r>100:
				n125=n125+1
			if r<150 and r>125:
				n150=n150+1
			if r<175 and r>150:
				n175=n175+1
			if r<200 and r>175:
				n200=n200+1 
			if r<225 and r>200:
				n225=n225+1
			if r<250 and r>225:
				n250=n250+1		
			if r<300 and r>250:
				n300=n300+1

	print (str(n50) +" --> "+ str(float(n50/1000))+"%")
	print (str(n75) +" -->"+ str(n75/1000)+"%") 
	print (str(n100)+" -->"+ str(n100/1000)+"%")
	print (str(n125)+" -->"+ str(n125/1000)+"%")
	print (str(n150)+" -->"+ str(n150/1000)+"%")
	print (str(n175)+" -->"+ str(n175/1000)+"%")
	print (str(n200)+" -->"+ str(n200/1000)+"%")
	print (str(n225)+" -->"+ str(n225/1000)+"%")
	print (str(n250)+" -->"+ str(n250/1000)+"%")
	print (str(n300)+" -->"+ str(n300/1000)+"%")
