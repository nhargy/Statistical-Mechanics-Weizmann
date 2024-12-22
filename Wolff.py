#######################
### WOLFF ALGORITHM ###
#######################

""" import libraries """
import numpy as np
from scipy.stats import moment
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
from labellines import labelLine, labelLines


###############
## functions ##
###############


def wolff_prob(b,j,h):
    return 1 - np.exp(-2*b*(j+h))

def toss(a):
    num = np.random.rand()
    if num < a:
        return True
    else:
        return False

def check(num1,num2,beta,j,h):
    if num1 == num2:
        if toss(wolff_prob(beta, j,-h*num2)) == True:
            return True
        return False
    else:
        return False

def lattice_gen(dim,random,arrow):
    lattice = []
    if random == True:
        for i in range(0,dim):
            row = []
            for j in range(0,dim):
                if toss(0.5) == True:
                    row.append(1)
                else:
                    row.append(-1)
            lattice.append(row)
    else:
         for i in range(0,dim):
            row = []
            for j in range(0,dim):
                row.append(arrow)
            lattice.append(row)
        
    return lattice


def cluster_gen(pos,lat,b,j,h,ant):
    cluster = [[pos[0],pos[1]]]
    N = len(lat[1])
    
    if ant == False:
    
        for i in cluster:
            x=i[0]; y=i[1]; curr=lat[x][y];
            if ([(x+1)%N,y] not in cluster):
                if check(curr, lat[(x+1)%N][y] , b,j,h) == True:
                    cluster.append([(x+1)%N,y])
            if ([(x-1)%N,y] not in cluster):
                if check(curr, lat[(x-1)%N][y] , b,j,h) == True:
                    cluster.append([(x-1)%N,y])
            if ([x,(y+1)%N] not in cluster):
                if check(curr, lat[x][(y+1)%N] , b,j,h) == True:
                    cluster.append([x,(y+1)%N])                
            if ([x,(y-1)%N] not in cluster):
                if check(curr, lat[x][(y-1)%N] , b,j,h) == True:
                    cluster.append([x,(y-1)%N])
                
    if ant == True:
        
        for i in cluster:
            x=i[0]; y=i[1]; curr=lat[x][y];
            if ([(x+1)%N,(y+1)%N] not in cluster):
                if check(curr, lat[(x+1)%N][(y+1)%N] , b,j,h) == True:
                    cluster.append([(x+1)%N,(y+1)%N])
            if ([(x-1)%N,(y+1)%N] not in cluster):
                if check(curr, lat[(x-1)%N][(y+1)%N] , b,j,h) == True:
                    cluster.append([(x-1)%N,(y+1)%N])
            if ([(x+1)%N,(y-1)%N] not in cluster):
                if check(curr, lat[(x+1)%N][(y-1)%N] , b,j,h) == True:
                    cluster.append([(x+1)%N,(y-1)%N])                
            if ([(x-1)%N,(y-1)%N] not in cluster):
                if check(curr, lat[(x-1)%N][(y-1)%N] , b,j,h) == True:
                    cluster.append([(x-1)%N,(y-1)%N])
    
    return cluster


def mag(Lattice,n):
    m = 0
    for i in range(0,n):
        for j in range(0,n):
            m+= Lattice[i][j]
    return m/(n*n)

def thermalise(Lattice,beta,J,h,MC,ant):
    size = len(Lattice[1])
    mean_mag = [mag(Lattice, size)]
    
    if ant == False:
    
        for n in range(0,MC):
            start = [np.random.randint(0,size),np.random.randint(0,size)]
            clust = cluster_gen(start, Lattice, beta, J,h,ant)
            for elm in clust:
                Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1)
            mean_mag.append(abs(mag(Lattice, size)))
            #plt.imshow(Lattice)
            #plt.show()
            
    if ant == True:
        
        for n in range(0,MC):
            start = [np.random.randint(0,size),np.random.randint(0,size)]
            start_1 = [(start[0]+1)%size,(start[1])%size]
            clust = cluster_gen(start, Lattice, beta, J,h,ant)
            clust_1 = cluster_gen(start_1, Lattice, beta, J,h,ant)
            for elm in clust:
                Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1)
            for elm in clust_1:
                Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1) 
            mean_mag.append(abs(mag(Lattice, size)))
            #plt.imshow(Lattice)
            #plt.show()
    
    
    return mean_mag, Lattice
    

def g(a,b):
    return 0.5 * (3 - (a/b**2))

############
## method ##
############


#simple phase transition



#init_lat = lattice_gen(50, True, 1)
beta_vals = np.linspace(0.3,0.7,35)

#beta_vals = [0.5]

sweeps = 40
back = 5
mags = []

for temp in beta_vals:
    init_lat = lattice_gen(30, True, 1)
    r = thermalise(init_lat, temp, 1, 0,sweeps,False)[1]
    avr = abs(np.mean(r))
    mags.append(avr)
 
grad = np.gradient(mags)
maxval = np.amax(grad)
index = np.where(grad == maxval)

b_val_max = beta_vals[index]

plt.plot(beta_vals,mags,color='blue', linewidth=2, marker='o',label=r'$m(\beta)$')   
plt.plot(beta_vals,grad,color='red', linewidth=2,linestyle='--',label=r'$dm/d \beta$')  
plt.axvline(b_val_max,color='black', linewidth=1, linestyle='dotted')
#labelLines(plt.gca().get_lines(), align=False, fontsize=11)
#plt.ylim(0.01)
plt.xlabel(r"$\beta$")
plt.ylabel(r"$m$")  
plt.title(r"Wolff Thermalisation: Ferromagnet ")
plt.grid("on")
plt.axvspan(b_val_max - 0.01, b_val_max + 0.01, color='green', alpha=0.5)
plt.legend()
plt.show()





#as a function of L



sweeps = 40
l_vals = [10,30,50]
beta_vals = np.linspace(0.35,0.6,35)

all_vals = []

for i in l_vals:
    g_vals = []
    for j in beta_vals:
        init_lat = lattice_gen(i, True, 1)
        R = thermalise(init_lat, j, 1, 0, sweeps, False)[1]
        av = abs(np.mean(R))
        g_vals.append(av)
    all_vals.append(g_vals)

grad1 = np.gradient(all_vals[0])
grad2 = np.gradient(all_vals[1])
grad3 = np.gradient(all_vals[2])

plt.plot(beta_vals,all_vals[0],label='L = {}'.format(str(l_vals[0])),marker='o',linewidth = 1.4,color='red')
plt.plot(beta_vals,all_vals[1],label='L = {}'.format(str(l_vals[1])),marker='o',linewidth = 1.4,color='blue')
plt.plot(beta_vals,all_vals[2],label='L = {}'.format(str(l_vals[2])),marker='o',linewidth = 1.4,color='black')

plt.plot(beta_vals,grad1,linewidth = 1,color='red',linestyle='--',)
plt.plot(beta_vals,grad2,linewidth = 1,color='blue',linestyle='--',)
plt.plot(beta_vals,grad3,linewidth = 1,color='black',linestyle='--',)


labelLines(plt.gca().get_lines(), align=False, fontsize=11)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$m$')
plt.title('Wolff Thermalisation: Ferromagnet')
plt.grid("on")
plt.show()



#anti-ferromagnet

beta_vals = np.linspace(0.3,0.7,25)

#beta_vals = [0.5]

L = 30
sweeps = 40
s = []

for temp in beta_vals:
    init_lat = lattice_gen(L, True, 1)
    a1 = thermalise(init_lat, temp, 1, 0,sweeps,False)[1]
    a2 = thermalise(init_lat, temp, 1, 0,sweeps,False)[1]
    mean_s = (abs(np.mean(a1)) - (-1)*abs(np.mean(a2)))/2
    s.append(mean_s)

grad = np.gradient(s)
maxval = np.amax(grad)
index = np.where(grad == maxval)

b_val_max = beta_vals[index]
 
plt.plot(beta_vals,s,color='green', linewidth=2, marker='o',label=r"$s$")
plt.plot(beta_vals,grad,color='red',linewidth=1.2,linestyle='--',label=r"$ds/d \beta$")    
plt.ylim(0.01)
plt.xlabel(r"$\beta$")
plt.ylabel(r"Staggered Magnetisation $s$")  
plt.title(r"Wolff Thermalisation: Anti-ferromagnet ")
plt.grid("on")
plt.axvspan(b_val_max - 0.01, b_val_max + 0.01, color='green', alpha=0.5)
plt.legend()
plt.show()


""" anti-ferromagnet with magnetic field """

B = -0.05

def mag_flip(beta,sig,H):
    return 1/(1+np.exp(-beta*sig*H))


def thermaliseH(Lattice,beta,J,h,MC,ant):
    size = len(Lattice[1])
    mean_mag = [mag(Lattice, size)]
    
    if ant == False:
    
        for n in range(0,MC):
            start = [np.random.randint(0,size),np.random.randint(0,size)]
            clust = cluster_gen(start, Lattice, beta, J,h,ant)
            if toss(mag_flip(beta,np.mean(clust),B))==False:
                for elm in clust:
                    Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1)
            mean_mag.append(abs(mag(Lattice, size)))
            #plt.imshow(Lattice)
            #plt.show()
            
    if ant == True:
        
        for n in range(0,MC):
            start = [np.random.randint(0,size),np.random.randint(0,size)]
            start_1 = [(start[0]+1)%size,(start[1])%size]
            clust = cluster_gen(start, Lattice, beta, J,h,ant)
            clust_1 = cluster_gen(start_1, Lattice, beta, J,h,ant)
            
            if toss(mag_flip(beta,np.mean(clust),B)==True):
                for elm in clust:
                    Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1)
            if toss(mag_flip(beta,np.mean(clust_1),B)==True):
                for elm in clust_1:
                    Lattice[elm[0]][elm[1]] = Lattice[elm[0]][elm[1]]*(-1) 
            mean_mag.append(abs(mag(Lattice, size)))
            #plt.imshow(Lattice)
            #plt.show()
    
    
    return mean_mag, Lattice
    

beta_vals = np.linspace(0.3,0.7,25)

L = 30
sweeps = 40
s = []

for temp in beta_vals:
    init_lat = lattice_gen(L, True, 1)
    init_lat2 = lattice_gen(L, True, 1)
    a1 = thermaliseH(init_lat, temp, 1, 0,sweeps,False)[1]
    a2 = thermaliseH(init_lat2, temp, 1, 0,sweeps,False)[1]
    mean_s = (abs(np.mean(a1)) - (-1)*abs(np.mean(a2)))/2
    s.append(mean_s)

grad = np.gradient(s)
maxval = np.amax(grad)
index = np.where(grad == maxval)

b_val_max = beta_vals[index]
 
plt.plot(beta_vals,s,color='green', linewidth=2, marker='o',label=r"$s$")
plt.plot(beta_vals,grad,color='red',linewidth=1.2,linestyle='--',label=r"$ds/d \beta$")    
plt.ylim(0.01)
plt.xlabel(r"$\beta$")
plt.ylabel(r"Staggered Magnetisation $s$")  
plt.title(r"Wolff Thermalisation: Anti-ferromagnet with B field ")
plt.grid("on")
plt.axvspan(b_val_max - 0.01, b_val_max + 0.01, color='green', alpha=0.5)
plt.legend()
plt.show()
