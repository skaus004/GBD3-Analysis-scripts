#pdb can be saved with 10 or 20 steps skipping; doesn't matter
#need to delete END from the pdb files using: grep -vwE "(END)" sourcefile > destinationfile
#command to run : python zaxis_percent.py filename.pdb number_of_frames output_file

import sys
import numpy as np
import linecache

outf = sys.argv[3]
info = open(outf,"w")

n=1
line_number=1
z_coord=np.array([])

pdb = sys.argv[1]                        #argument for command line
frames = sys.argv[2]
info.write(pdb+" with frames "+frames+'\n')
frames = int(frames)
with open(pdb) as f1:
    next(f1)                             #skip first line

    current_line = 1
    for line in f1:
        if current_line == line_number:  #matching lines with the required line_number
            z=(line[47:50])
            x=(line[30:34])
#           x=(line[38:42])
            zplusx=(str(z)+","+str(x))
            z_coord=np.append(z_coord,zplusx) #listing all x coordinates
            n=n+84
            line_number = n              #changing required line_number from 1 to 85 to read next molecule
        current_line += 1                #reading file line by line to match with the line number in the next loop
z_coord = np.reshape(z_coord,(frames,-1))   #shaping x_coord array having columns as frames,i.e. 22 is number of frames
z_coordt=z_coord.transpose()
print (z_coordt)
z0=0
z20=0
z40=0
z60=0
z80=0
z100=0
z120=0
z140=0
z160=0
z180=0
z200=0
Az0=np.array([])
Az20=np.array([])
Az40=np.array([])
Az60=np.array([])
Az80=np.array([])
Az100=np.array([])
Az120=np.array([])
Az140=np.array([])
Az160=np.array([])
Az180=np.array([])
Az200=np.array([])

for row in z_coord:
    for cell in row:
        z=int(cell[:3])
        x=int(cell[4:])
        x2=x*x
        if x2<22500:  
            if z<0:
                z0=z0+1
            if (z<20 and z>0):
                z20=z20+1
            if z<40 and z>20: 
                z40=z40+1
            if z<60 and z>40:
                z60=z60+1
            if z<80 and z>60:
                z80=z80+1
            if z<100 and z>80:
                z100=z100+1
            if z<120 and z>100:
                z120=z120+1
            if z<140 and z>120:
                z140=z140+1
            if z<160 and z>140:
                z160=z160+1
            if z<180 and z>160:
                z180=z180+1
            if z<200 and z>180:
                z200=z200+1


print ("g0="+str(z0/((100*0.068182)*frames)))
print ("g20="+str(z20/((100*0.068182)*frames)))
print ("g40="+str(z40/((100*0.068182)*frames)))
print ("g60="+str(z60/((100*0.068182)*frames)))
print ("g80="+str(z80/((100*0.068182)*frames)))
print ("g100="+str(z100/((100*0.068182)*frames)))
print ("g120="+str(z120/((100*0.068182)*frames)))
print ("g140="+str(z140/((100*0.068182)*frames)))
print ("g160="+str(z160/((100*0.068182)*frames)))
print ("g180="+str(z180/((100*0.068182)*frames)))
print ("g200="+str(z200/((100*0.068182)*frames)))

'''
    totalz=z0+z20+z40+z60+z80+z100+z120+z140+z160+z180+z200
    pz0=(z0*100)/totalz
    pz20=(z20*100)/totalz
    pz40=(z40*100)/totalz
    pz60=(z60*100)/totalz
    pz80=(z80*100)/totalz
    pz100=(z100*100)/totalz
    pz120=(z120*100)/totalz
    pz140=(z140*100)/totalz
    pz160=(z160*100)/totalz
    pz180=(z180*100)/totalz
    pz200=(z200*100)/totalz

    Az0=np.append(Az0,pz0)
    Az20=np.append(Az20,pz20)
    Az40=np.append(Az40,pz40)
    Az60=np.append(Az60,pz60)
    Az80=np.append(Az80,pz80)
    Az100=np.append(Az100,pz100)
    Az120=np.append(Az120,pz120)
    Az140=np.append(Az140,pz140)
    Az160=np.append(Az160,pz160)
    Az180=np.append(Az180,pz180)
    Az200=np.append(Az200,pz200)

print("z<0   ---> mean="+str(np.mean(Az0))+", sd = "+str(((np.std(Az0))*100)/(np.mean(Az0))))
print("z<20  ---> mean="+str(np.mean(Az20))+", sd = "+str(((np.std(Az20))*100)/(np.mean(Az20))))
print("z<40  ---> mean="+str(np.mean(Az40))+", sd = "+str(((np.std(Az40))*100)/(np.mean(Az40))))
print("z<60  ---> mean="+str(np.mean(Az60))+", sd = "+str(((np.std(Az60))*100)/(np.mean(Az60))))
print("z<80  ---> mean="+str(np.mean(Az80))+", sd = "+str(((np.std(Az80))*100)/(np.mean(Az80))))
print("z<100 ---> mean="+str(np.mean(Az100))+", sd = "+str(((np.std(Az100))*100)/(np.mean(Az100))))
print("z<120 ---> mean="+str(np.mean(Az120))+", sd = "+str(((np.std(Az120))*100)/(np.mean(Az120))))
print("z<140 ---> mean="+str(np.mean(Az140))+", sd = "+str(((np.std(Az140))*100)/(np.mean(Az140))))
print("z<160 ---> mean="+str(np.mean(Az160))+", sd = "+str(((np.std(Az160))*100)/(np.mean(Az160))))
print("z<180 ---> mean="+str(np.mean(Az180))+", sd = "+str(((np.std(Az180))*100)/(np.mean(Az180))))
print("z<200 ---> mean="+str(np.mean(Az200))+", sd = "+str(((np.std(Az200))*100)/(np.mean(Az200))))
'''


