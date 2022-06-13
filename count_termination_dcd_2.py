import dcd as dcd
import numpy as np
import sys
import math

if len(sys.argv)<2:
        print("Usage --> python count_termination_dcd.py <input dcd> <index of center of mass of ligand> <no of replicas> <dinstance between initial and binding sites in Angstrom> <x-coordinate of direct binding box center> <y...> <z...> <direct binding box dimension x> <y...> <z...> <x-coordinate of dwelling sphere> <y...> <z...> <dwell sphere radius>")


if len(sys.argv)==15:
	input_dcd = sys.argv[1]
	com_index =int(sys.argv[2])
	n_replicas = int(sys.argv[3])
	binding_jump = int(float(sys.argv[4]))
	c1 = int(float(sys.argv[5]))
        c2 = int(float(sys.argv[6]))
        c3 = int(float(sys.argv[7]))
        direct_x_cut_off=  int(float(sys.argv[8]))
	direct_y_cut_off= int(float(sys.argv[9]))
	direct_z_cut_off= int(float(sys.argv[10]))
	d1 = int(float(sys.argv[11]))
	d2 = int(float(sys.argv[12]))
	d3 = int(float(sys.argv[13]))
	dwell_cut_off = int(float(sys.argv[14]))

	dind_x1= (direct_x_cut_off/2)+c1
	dind_x2= (direct_x_cut_off/2)-c1
	dind_y1= (direct_y_cut_off/2)+c2
	dind_y2= (direct_y_cut_off/2)-c2
	dind_z1= (direct_z_cut_off/2)+c3
	dind_z2= (direct_z_cut_off/2)-c3


	print ("Center of mass index     = "+str(com_index))
	print ("no. of replicas          = "+ str(n_replicas))
	cc,co=dcd.read_dcd_data_ref(input_dcd)  #co are the cordinates in numpy array [frames, atoms, xyz]==[axis0, axis1, axis2]
	n_atm_per_res=int(int(np.size(co, 1))/n_replicas)

	print ("no. of atoms per residue = "+ str(n_atm_per_res))
	print ("no. of frames            = "+str(np.size(co,0)))
	print ("     ")
	print (np.shape(co))

	co=co[:, com_index:-1:n_atm_per_res, :]  # [frames, com_atom, xyz] == [axis0, axis1, axis2]
	print (np.shape(co))

	####calculating r2-r1 for each step in simulation

	co_1 = np.delete(co,0, axis=0) #deleting first row in co means deleting first frame
	co_2 = np.delete(co,-1, axis=0) #deleting last row in co means deleting last frame

	co1_co2=np.subtract(co_1,co_2)
	co1_co2_sq=np.power(co1_co2,2)
	r1_r2_sq=np.sum(co1_co2_sq, axis=2)
	r1_r2=np.power(r1_r2_sq,0.5)
	r1_r2=r1_r2.astype('int') 
	r1_r2_sq=r1_r2_sq.astype('int')

	####calculating r for each step in simulation

	co_sq=np.power(co,2)   #x^2, y^2, z^2
	r2=np.sum(co_sq, axis=2)  #r2=(x*x + y*y + z*z)
	r=np.power(r2,0.5)  # data of r in form of [frames, com_atom]
	r=r.astype('int')
	all_co_r=np.dstack((co, r))  # [frames, com_atom, x y z r] == [axis0, axis1, axis2]
	all_co_r1_r2_sq=np.dstack((co_1,r1_r2_sq))
 
	binding_jump_sq=binding_jump*binding_jump
	N_termination = r1_r2_sq[r1_r2_sq>binding_jump_sq]
	total_binds = int(np.size(N_termination))+1

	####direct bindings
	b1=0
	b2=0
	b3=0
	n_slicing=np.array([])
	m_slicing=np.array([])
	n_frames = int(np.size(co,0))-1
	for m in range(0,n_replicas):
		for n in range(0,n_frames):
			if all_co_r1_r2_sq[n][m][3]>binding_jump_sq:
				n_slicing=np.append(n_slicing,n)
				m_slicing=np.append(m_slicing,m)

	slicing=np.vstack((n_slicing,m_slicing))   # 2d array having info about which replicas binds at which frame [2 rows: 1st row is frame number and 2nd row is replica number ] 

	for p in range(0,np.size(slicing,1)):
		if p==0:
                        f2=int(slicing[0][p])
                        rep_no=int(slicing[1][p])
                        temp=all_co_r1_r2_sq[:f2,rep_no,:]   # 2d array having axis0 as frames and axis1 as x y z r^2 for each replica simulation before termination
                        for q in range(0,np.size(temp,0)):
                                x=int(temp[q][0])
                                y=int(temp[q][1])
                                z=int(temp[q][2])
                                if (dind_x1<x or x>dind_x2 or dind_y1<y or y>dind_y2 or dind_z1<z or z>dind_z2):                              
                               #direct_r =  math.sqrt((x-c1)*(x-c1) + (y-c2)*(y-c2) + (z-c3)*(z-c3))
                               # if direct_r > direct_cut_off:
                                        b1=b1+1
                                        break

		if slicing[1][p]==slicing[1][p-1]:
			f1=int(slicing[0][p-1])
			f2=int(slicing[0][p])
			rep_no=int(slicing[1][p])
			temp=all_co_r1_r2_sq[f1:f2,rep_no,:]
                        for q in range(0,np.size(temp,0)):
                                x=int(temp[q][0])
                                y=int(temp[q][1])
                                z=int(temp[q][2])
                                if (dind_x1<x or x>dind_x2 or dind_y1<y or y>dind_y2 or dind_z1<z or z>dind_z2):
                                # direct_r =  math.sqrt((x-c1)*(x-c1) + (y-c2)*(y-c2) + (z-c3)*(z-c3))
                               # if direct_r > direct_cut_off:
                                        b2=b2+1
                                        break
			
		if slicing[1][p]!=slicing[1][p-1]:	
			f2=int(slicing[0][p])
			rep_no=int(slicing[1][p])
			temp=all_co_r1_r2_sq[:f2,rep_no,:]
			
			for q in range(0,np.size(temp,0)):
				x=int(temp[q][0])
				y=int(temp[q][1])
				z=int(temp[q][2])
                                if (dind_x1<x or x>dind_x2 or dind_y1<y or y>dind_y2 or dind_z1<z or z>dind_z2):
#				direct_r =  math.sqrt((x-c1)*(x-c1) + (y-c2)*(y-c2) + (z-c3)*(z-c3))  
#				if direct_r > direct_cut_off:
					b3=b3+1
					break
	print ("     ")						
	print ("indirect bindings ---> "+ str(b1+b2+b3))
	print ("direct bindings   ---> "+ str(total_binds-b1-b2-b3))			
	print ("total bindings    ---> "+ str(total_binds))	
	
	####dwelling time
        s=0
        for t in range(0,np.size(slicing,1)):
                if p==0:
                        f2=int(slicing[0][p])
                        rep_no=int(slicing[1][p])
                        temp=all_co_r1_r2_sq[:f2,rep_no,:]
                        for q in range(-1*int(np.size(temp,0))+1,1):
                                q=-1*q
                                x=int(temp[q][0])
                                y=int(temp[q][1])
                                z=int(temp[q][2])
                                dwell_r =  math.sqrt((x-d1)*(x-d1) + (y-d2)*(y-d2) + (z-d3)*(z-d3))
                                if dwell_r > dwell_cut_off:
                                        s=s+int(np.size(temp,0))-q
                                        break

                if slicing[1][p]==slicing[1][p-1]:
                        f1=int(slicing[0][p-1])
                        f2=int(slicing[0][p])
                        rep_no=int(slicing[1][p])
                        temp=all_co_r1_r2_sq[f1:f2,rep_no,:]
                        for q in range(-1*int(np.size(temp,0))+1,1):
                                q=-1*q
                                x=int(temp[q][0])
                                y=int(temp[q][1])
                                z=int(temp[q][2])
                                dwell_r =  math.sqrt((x-d1)*(x-d1) + (y-d2)*(y-d2) + (z-d3)*(z-d3))
                                if dwell_r > dwell_cut_off:
                                        s=s+int(np.size(temp,0))-q
                                        break

                if slicing[1][p]!=slicing[1][p-1]:
                        f2=int(slicing[0][p])
                        rep_no=int(slicing[1][p])
                        temp=all_co_r1_r2_sq[:f2,rep_no,:]

                        for q in range(-1*int(np.size(temp,0))+1,1):
				q=-1*q
                                x=int(temp[q][0])
                                y=int(temp[q][1])
                                z=int(temp[q][2])
                                dwell_r =  math.sqrt((x-d1)*(x-d1) + (y-d2)*(y-d2) + (z-d3)*(z-d3))
                                if dwell_r > dwell_cut_off:
                                        s=s+int(np.size(temp,0))-q
                                        break
	print ("Avg dwell time    ---> "+str(s/total_binds)+"frames * <writeTraj> * <timestep>")

