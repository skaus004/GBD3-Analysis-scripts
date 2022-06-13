#Usage--> python n_binds.py <logfile>

import sys
import numpy as np
import matplotlib.pylab as plt
log_file = sys.argv[1]
f = open(log_file)

last_n_elements = 100
n=0

B_T = np.array([])
Number_of_bindings = np.array([])

for line in f :
    if line[0:2] == "#1" :
        n=n+1
        time=float((line.split() [6]).replace('t=',''))
        #print ("time = "+str(time))
        B_T = np.append(B_T, time) #B_T = [t1 , t2 , t3 , ...]
        Number_of_bindings = np.append(Number_of_bindings,n) # Number_of_bindings = [1,2,3,...]
        if n == 600 : #giving data upto 600 bindings
       	    break
    elif line[2:6] == "Step" :
        step = line.split()[2]
        full_step = line[0:30]


B_T_progression = np.cumsum(B_T) #B_T_progression = [t1 , t1+t2 , t1 + t2 + t3 ,...]
T_average = np.divide(B_T_progression,Number_of_bindings) # T_average = [t1/1 , (t1+t2)/2 , (t1+t2+t3)/3 ,...]
last_100_T_average = T_average[-100:] # last_100_T_average = [array of last 100 elements of T_average] for standard deviation calculation
sd_t_avg = np.std(last_100_T_average)

print ("relative sd of t_avg for last 100 bindings = " + str(((sd_t_avg*100)/last_100_T_average.mean())) + " %")
print ("k_ON = " + str(1000000/(1.1*(B_T.mean()))) + " e^10 M^-1 s^-1")
print ("last step = "+str(full_step))
print ("running time = "+ str(abs(float(step))/(1000000*20))+ " microsec")
print ("Bound ligands = "+str(n))
print ("Average binding time = " + str(B_T.mean()/1000000) + "microsec")
print ("maximum binding time = " + str(B_T.max()) + "ps")
plt.plot(T_average)
plt.ylabel('Binding time average (ps)')
plt.show()

