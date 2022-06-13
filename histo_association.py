import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)==3:
    log_file = sys.argv[1]
    log_file2 =sys.argv[2]
    f = open(log_file)
    f2 = open(log_file2)
    n=0
    B_T = np.array([])
    Number_of_bindings = np.array([])
    wrd='satisfied at t='
    for line in f :
        if wrd in line:
            all_time= line.split(wrd,1)[1]
            time=float(all_time.split()[0])/1000000
            n=n+1
            B_T = np.append(B_T, time) #B_T = [t1 , t2 , t3 , ...]
            Number_of_bindings = np.append(Number_of_bindings,n) # Number_of_bindings = [1,2,3,...]

    for line2 in f2 :
        if wrd in line2:
            all_time=line2.split(wrd,1)[1]
            time=float(all_time.split()[0])/1000000
            n=n+1
            B_T=np.append(B_T, time)
            Number_of_bindings = np.append(Number_of_bindings,n)
            if n==600:
                break

if len(sys.argv)==2:
    log_file = sys.argv[1]
    f=open(log_file)
    n=0
    B_T =np.array([])
    Number_of_bindings = np.array([])
    wrd='satisfied at t='
    for line in f:
        if wrd in line:
            all_time=line.split(wrd,1)[1]
            time=float(all_time.split()[0])/1000000
            n=n+1
            B_T = np.append(B_T, time)
            Number_of_bindings = np.append(Number_of_bindings,n)
            if n==600:
                break
print ("Bound ligands = "+str(n))
print (str(np.shape(B_T)))
plt.hist(B_T, bins=100)
plt.xlabel('Binding time (microsec)')
plt.ylabel('Number of bindings')
plt.savefig('figure.png')

