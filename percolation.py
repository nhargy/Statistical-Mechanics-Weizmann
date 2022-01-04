#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import time
get_ipython().run_line_magic('matplotlib', 'inline')


# In[398]:


### functions ###

def prob(p):
    if np.random.rand() <= p:
        return True
    else:
        return False

def lattice_gen(N,p):
    lattice = []
    count_a = 0  # keeps track of even/odd row
    for i in range(0,N):
        row = []
        if np.mod(count_a,2) == 0: # an even row is a row with no zeros
            count_b = 0            # keeps track of even/odd column
            for j in range(0,N):
                if np.mod(count_b,2) == 0:
                    row.append(1)
                else:
                    if prob(p) == True:
                        row.append(1)
                    else:
                        row.append(5)
                count_b+=1
        else:                     # an odd row is a row with zeros every other index
            count_b = 0           # keeps track of even/odd column
            for j in range(0,N):
                if np.mod(count_b,2) == 0:
                    if prob(p) == True:
                        row.append(1)
                    else:
                        row.append(5)
                else:
                    row.append(0)
                count_b+=1
        lattice.append(row)
        count_a+=1
    return lattice

def DFS(lattice):
    N = len(lattice[0])
    path = []
    great_success = 0
    
    for i in range(0,N):
        if lattice[0][i] == 1:
            path.append([0,i])
    
    visited = []
    for pos in path:
        if pos not in visited:
            x = pos[0]; y = pos[1]; index = path.index(pos);
            
            if lattice[N-1][1]==0:
                check = N-2
            else:
                check= N-1
            
            if (x == check):  # or (x == N -2) :
                great_success = 1
                break
            visited.append(pos)

            if x != 0:
                if lattice[x-1][y] == 1:
                    path.insert(index+1,[x-1,y])
            if lattice[x][(y+1)%N] == 1:
                path.insert(index+1,[x,(y+1)%N])
            if lattice[x][(y-1)%N] == 1:
                path.insert(index+1,[x,(y-1)%N])
            if lattice[x+1][y] == 1:
                path.insert(index+1,[x+1,y])
    
    return great_success

def decimate(lattice):
    N = len(lattice[0])
    new_lattice = []
    left_cords = [[0,1],[0,3],[0,4],[1,4],[3,4]]
    
    row_count = 0
    
    newr_count = 0
    
    while row_count < N:
        column_count = 0
        new_lattice.append([]); new_lattice.append([]);
        while column_count < N:
            new_box = []
            
            for i in range(row_count,row_count+5):
                roo = []
                for j in range(column_count,column_count+5):
                    roo.append(lattice[i%N][j%N])
                new_box.append(roo)
            left_box = []
            
            for i in range(0,5):
                left_row = []
                for j in range(0,5):
                    if [i,j] in left_cords:
                        left_row.append(0)
                    else:
                        left_row.append(new_box[i][j])
                        
                left_box.append(left_row)
            left_bottom = DFS(left_box)
            left_left = DFS(np.transpose(left_box))

            
            dec_box = [[1,0],[1,1]]
            
            if left_left == 0:
                dec_box[1][1] = 5
            if left_bottom == 0:
                dec_box[0][0] = 5

            for i in range(0,2):
                for j in range(0,2):
                    new_lattice[newr_count + i].append(dec_box[i][j])
            
            column_count+=4
        newr_count+=2
        row_count+=4
     
    return new_lattice


# In[402]:


n = 259
resolution = 20
p_vals = np.linspace(0.01,1,resolution)
sweeps = 100
steps = sweeps*resolution

P_vals = []
load = 0
for p in p_vals:
    mean = 0
    for j in range(0,sweeps):
        print(str(r'Loading (p={}): ').format(str(round(p,2)))+str(round((load/steps)*100,1))+'%',end='\r')
        box = lattice_gen(n,p)
        mean+=DFS(box)
        load+=1
    P_vals.append(mean/sweeps)

gradient = np.gradient(P_vals)    

plt.plot(p_vals, P_vals,color='blue',marker='.',linewidth = 2.3,label='P')
plt.plot(p_vals, gradient,linestyle='--',color='red',label='dP/dp')

plt.grid('on')
plt.xlabel(r'$p$')
plt.ylabel('Probability of Path')
plt.title('Probability of Percolation')
plt.legend()
plt.show()


# In[396]:


p_vals = [0.4,0.5,0.6]
N = 80
n = 50

steps = len(p_vals)*n
load = 0
for p in p_vals:
    s = 0
    s1 = 0
    s2 = 0
    s3 = 0
    for i in range(0,n):
        print(str(r'Loading (p={}): ').format(str(round(p,2)))+str(round((load/steps)*100,1))+'%',end='\r')
        box = lattice_gen(300,p)
        box1 = decimate(box)
        box2 = decimate(box1)
        box3 = decimate(box2)
        
        s+=DFS(box)
        s1+=DFS(box1)
        s2+=DFS(box2)
        s3+=DFS(box3)
        
        load+=1
    
    S = s/n
    S1 = s1/n
    S2 = s2/n
    S3 = s3/n
    
    P = [S,S1,S2,S3]
    plt.plot(P,label='p={}'.format(str(p)),marker='o')
plt.xlabel("Decimation Steps")
plt.ylabel("Probability of Path")
plt.grid("on")
plt.legend()
plt.show()
