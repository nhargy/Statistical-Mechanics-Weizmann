#!/usr/bin/env python
# coding: utf-8

# In[216]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns; sns.set_theme()


# In[197]:


########################
### DEFINE FUNCTIONS ###
########################

def disorder_gen(mean,SD,N):       ### CREATES DISORDER MATRICES H AND J
    mat = np.random.normal(mean,SD,size=N)
    try:
        for i in range(0,N[0]):
            mat[i][i] = 0
    except TypeError:
        pass
    return mat.round(5)

def lattice_gen(N):                ### GENERATES 1D RING OF SPINS
    lattice = []
    for i in range(0,N):
        a = np.random.rand()
        if a < 0.5:
            lattice.append(-1)
        else:
            lattice.append(1)
    return lattice

def ham(J,h,N,s1,s2,model):         ### HAMILTONIAN
    if model=='RFIM':
        return - ((J/N)*s1*s2) - (h*s1)
    elif model=='SK':
        return - ((J/N)*s1*s2)
    

def thermalise(lattice,J,h,beta,model,sweeps):     ### METROPOLIS
    N = len(lattice)
    J_const = np.mean(J)
    
    for s in range(0,sweeps):
        pos = np.random.randint(0,N)
        spin = lattice[pos]
        f_spin = spin*(-1)
        
        if model=='RFIM':
            energy0 = 0
            energy1 = - ham(J_const,h[pos],N,f_spin,lattice[pos],model)
            for i in range(0,N):
                a0 = ham(J_const,h[pos],N,spin,lattice[i],model)
                a1 = ham(J_const,h[pos],N,f_spin,lattice[i],model)
                energy0+=a0
                energy1+=a1
                
            if energy1 < energy0:
                lattice[pos] = f_spin
            else:
                diff = energy1-energy0
                if np.random.rand() < np.exp(-beta*diff):
                    lattice[pos] = f_spin
    
        if model=='SK':
            energy0 = 0
            energy1 = 0
            for i in range(0,N):
                a0 = ham(J[pos][i],0,N,spin,lattice[i],model)
                a1 = ham(J[pos][i],0,N,f_spin,lattice[i],model)
                energy0+=a0
                energy1+=a1
                
            if energy1 < energy0:
                lattice[pos] = f_spin
            else:
                diff = energy1-energy0
                if np.random.rand() < np.exp(-beta*diff):
                    lattice[pos] = f_spin           
    
    return lattice

def Q_par(array):
    Q = []
    spins = len(array[0])
    replicas = len(array)
    
    for a in range(0,replicas):
        row = []
        for b in range(0,replicas):
            q_ab = 0
            for i in range(0,spins):
                val = array[a][i]*array[b][i]
                q_ab+=val
            row.append((1/spins)*q_ab)
        Q.append(row)
    return Q
        


# In[328]:


###########################################################
### CREATE AND THERMALISE REPLICAS FOR BOTH RFIM AND SK ###
###########################################################

# CONSTANTS
N = 500
beta = 150
mean = 0
sweeps = 1000 

# DISORDER (FROZEN)
J = disorder_gen(mean,N,(N,N))
h = disorder_gen(mean,1,N)

runs = 400
array_RFIM = []
load = 0

for r in range(0,runs):
    print(str(r'Loading RFIM model (replica={}): ').format(str(r))+str(round((load/runs)*100,1))+'%',end='\r')
    box = lattice_gen(N)
    g = thermalise(box,J,h,beta,'RFIM',sweeps)
    array_RFIM.append(g)
    load+=1
    
array_SK = []
load = 0

for r in range(0,runs):
    print(str(r'Loading SK model (replica={}): ').format(str(r))+str(round((load/runs)*100,1))+'%',end='\r')
    box = lattice_gen(N)
    g = thermalise(box,J,h,beta,'SK',sweeps)
    array_SK.append(g)
    load+=1


# In[329]:


##########################################
### CALCULATE ORDER-PARAMETER MATRICES ###
##########################################

Q_RFIM = Q_par(array_RFIM)
Q_SK = Q_par(array_SK)


# In[330]:


##########################
### PLOT HEATMAPS OF Q ###
##########################

ax = plt.axes()
sns.heatmap(Q_RFIM,vmin=-1, vmax=1,ax=ax, xticklabels=50, yticklabels=False)

ax.set_title(r'$Q_{\alpha \beta}$  Random Field Ising')
plt.show()

ax = plt.axes()
sns.heatmap(Q_SK,vmin=-1, vmax=1,ax=ax, xticklabels=50, yticklabels=False)

ax.set_title(r'$Q_{\alpha \beta}$  Sherington-Kirkpatrick')
plt.show()


# In[331]:


Q_RFIM_LIST = []
ni = len(Q_RFIM[0])
for i in range(0,ni):
    for j in range(0,ni):
        Q_RFIM_LIST.append(Q_RFIM[i][j])
        
Q_SK_LIST = []
ni = len(Q_SK[0])
for i in range(0,ni):
    for j in range(0,ni):
        Q_SK_LIST.append(Q_SK[i][j])

sns.histplot(Q_RFIM_LIST,binwidth=0.008,kde=True)
sns.histplot(Q_SK_LIST,binwidth=0.008,kde=True)
plt.xlabel(r'$Q_{\alpha \beta}$')
plt.title('Histogram: RFIM vs SK')
plt.xlim(-1,1)
plt.annotate("RFIM", (0.27, 6000))
plt.annotate("SK", (-0.22, 6000))


# In[ ]:




