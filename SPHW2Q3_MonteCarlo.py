###########################
### MONTE CARLO METHODS ###
### By Nadav Hargittai  ###
###########################

""" Import libraries """
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines

""" Define cosntants """
J = 1
B = 1
N = 512

""" Define plotting constants """

W = 2.5
F = 15

#####################
# 111 Functions 111 #
#####################

""" 1) Hamiltonian function: gives either one spin flip energy difference or
total energy of input lattice """
def H(couple,magnet,temp,spins,position,local):
    Hml = lambda a,b : -J*a*b - B*a
    if local==True:
        H_current = Hml(spins[(position-1)%N],spins[position%N]) + Hml(spins[(position)%N],spins[(position+1)%N])
        H_new = Hml(spins[(position-1)%N],(-1)*spins[position%N]) + Hml((-1)*spins[(position)%N],spins[(position+1)%N])
        
        delta = H_new - H_current
        
        if (H_new - H_current) < 0:
            return True, delta
        else:
            return False, delta
        
    else:
        H_sum = 0
        H_sum2 = 0
        for i in range(0,N):
            z = Hml(spins[i%N],spins[(i+1)%N])
            H_sum += z
            H_sum2 += z**2
        
        return H_sum/len(spins), H_sum2/len(spins)

""" 2) Mean magnetisation function """
def MAG(spins):
    size = len(spins)
    m = 0
    for i in spins:
            m+=i
    return m/size

""" 3) Lattice creation function """

def lattice(size,random):
    lat = []
    if random == True:
        for i in range(0,size+1):
            a = np.random.rand()
            if a > 0.5:
                lat.append(1)
            elif a < 0.5:
                lat.append(-1)
            else:
                lat.append(1)
        return lat
    if random == False:
        for i in range(0,size+1):
            lat.append(1)
        return lat

""" 4) Thermalisation function: takes a lattice as input and produces the
thermalised state. Returns the magentisation m as a function of MC sweeps. """
def thermalise(spins,MC,sweeps,beta,method):
    if method == 'M':
        m = [MAG(spins)]
        for i in range(0,MC):
            for j in range(0,sweeps):
                flip = round(np.random.rand()*N) % N
                energy = H(J,B,beta,spins,flip,True)
                if energy[0] == True:
                    spins[flip%N] = spins[flip%N]*(-1)
                elif energy[0] == False:
                    boltz = np.exp(-beta*energy[1])
                    flip1 = np.random.rand()
                    if (flip1 < boltz):
                        spins[flip%N] = spins[flip%N]*(-1)
            m.append(MAG(spins))
        
    if method == 'G':
            m=[MAG(spins)]
            for i in range(0,MC):
                for j in range(0,sweeps):
                    flip = round(np.random.rand()*N) % N
                    energy = H(J,B,beta,spins,flip,True)
    
                    boltz = 1/(1+np.exp(-beta*energy[1]))
                    flip1 = np.random.rand()
                    if (flip1 > boltz):
                        spins[flip%N] = spins[flip%N]*(-1)
                m.append(MAG(spins))
    
    return m,spins

""" 5) Correlation as a function of distance """
def corr_as_k(spins):
    run = int(round(N/2))
    k_vals = np.arange(1,run-1)
    corr_vals = []
    for i in k_vals:
        tot = 0
        for j in range(0,N):
            a = spins[int(j)%N]*spins[int(j+i)%N]
            tot+=a
        corr_vals.append((tot/N))
    return corr_vals



##################
# 222 Method 222 #
##################

MC_sweeps = 35

""" 1) Thermalise a random lattice for three representative betas """
a1 = thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[0]
a2 = thermalise(lattice(N,True), MC_sweeps, N, 0.5, 'M')[0]
a3 = thermalise(lattice(N,True), MC_sweeps, N, 0.2, 'M')[0]

plt.plot(a1,label=r'$\beta = 1$',color = 'red',linewidth=W); 
plt.plot(a2,label=r'$\beta = 0.5$',color = 'blue',linewidth=W); 
plt.plot(a3,label=r'$\beta = 0.2$',color = 'black',linewidth=W)
labelLines(plt.gca().get_lines(), align=False, fontsize=F)
plt.title('Metropolis Thermalisation: Random Lattice Initial Conditions')
plt.xlabel('MC Sweeps')
plt.ylabel('Mean Magnetisation')
plt.grid('on')
plt.show()


""" 2) Thermalise an all-up lattice for three representative betas """
b1 = thermalise(lattice(N,False), MC_sweeps, N, 1, 'M')[0]
b2 = thermalise(lattice(N,False), MC_sweeps, N, 0.5, 'M')[0]
b3 = thermalise(lattice(N,False), MC_sweeps, N, 0.2, 'M')[0]

plt.plot(b1,label=r'$\beta = 1$',color = 'red',linewidth=W); 
plt.plot(b2,label=r'$\beta = 0.5$',color = 'blue',linewidth=W); 
plt.plot(b3,label=r'$\beta = 0.2$',color = 'black',linewidth=W)
labelLines(plt.gca().get_lines(), align=False, fontsize=F)
plt.title('Metropolis Thermalisation: All-up Lattice Initial Conditions')
plt.xlabel('MC Sweeps')
plt.ylabel('Mean Magnetisation')
plt.grid('on')
plt.show()

""" 3) Thermalise both Glabuer and Metropolis, and compare """
c1 = thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[0]
c2 = thermalise(lattice(N,True), MC_sweeps, N, 1, 'G')[0]
plt.plot(c1,label='Metropolis',color = 'red',linewidth=W-1); 
plt.plot(c2,label='Glauber',color = 'green',linewidth=W-1); 
plt.legend()
plt.title('Metropolis vs Glauber')
plt.xlabel('MC Sweeps')
plt.ylabel('Mean Magnetisation')
plt.grid('on')
plt.show()



######################
# 333 Properties 333 #
######################


""" 1) Calculate the internal energy of a thermalised system (Metropolis) 
many times to get a good mean"""

e_vals = []
for i in range(0,100):
    run = thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[1]
    U = H(J,B,1,run,1,False)[0]
    e_vals.append(U)

plt.hist(e_vals,color='red')
plt.title('Mean Energy Frequency')
plt.xlabel('Energy per spin')
plt.ylabel('Frequency')
plt.show()

mean_energy = np.mean(e_vals)
print('The mean energy per spin is '+str(mean_energy))


""" 2) Calculate the mean of the energy squared """

e_vals2 = []
for i in range(0,100):
    run = thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[1]
    U = H(J,B,1,run,1,False)[1]
    e_vals2.append(U)

mean_energy2 = np.mean(e_vals2)
print('The squared energy mean: ' +str(mean_energy2))

cv = (mean_energy2 - mean_energy**2)
print('The specific heat is: '+str(cv))


""" 3) Calculate the correlation """

d1 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[1])
d2 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 0.5, 'M')[1])
d3 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 0.2, 'M')[1])
plt.plot(d1,label=r'$\beta = 1$',color='red')
plt.plot(d2,label=r'$\beta = 0.5$',color='blue')
plt.plot(d3,label=r'$\beta = 0.2$',color='black')
plt.title('Correlation as a function of distance: B field on')
plt.xlabel('Distance k')
plt.ylabel('Correlation')
plt.grid('on')
plt.legend()
plt.show()


B = 0
e1 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 1, 'M')[1])
e2 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 0.5, 'M')[1])
e3 = corr_as_k(thermalise(lattice(N,True), MC_sweeps, N, 0.2, 'M')[1])
plt.plot(e1,label=r'$\beta = 1$',color='red')
plt.plot(e2,label=r'$\beta = 0.5$',color='blue')
plt.plot(e3, label=r'$\beta = 0.2$',color='black')
plt.title('Correlation as a function of distance: B field off')
plt.xlabel('Distance k')
plt.ylabel('Correlation')
plt.grid('on')
plt.legend()
plt.show()
