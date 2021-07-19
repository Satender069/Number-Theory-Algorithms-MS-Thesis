#!/usr/bin/env python
# coding: utf-8

# ###  n : number to be tested,  k : number of times test runs before returning "Probable Prime"

# In[45]:


import random 
def fermat(n,k):
    for i in range(k):
        a = random.randint(2,n)
        if pow(a,n,n)!= pow(a,1,n):
            return "Definite Composite"
    return "Probable Prime"


# In[46]:


fermat(163751584889877593680759288789633023806794889209326665885982508719719335997517,5)


# In[47]:


fermat(232354763543765732524573257,5)


# In[48]:


fermat(8489327852164434890945021112361,15)


# #### Time analysis

# Here we show how time of the execution of above program varies. We check running time of the program for a random $n$-$digit$ number where $2\leq n\leq B$. $B$ is given as an input.

# In[49]:


import matplotlib.pyplot as plt
import time
def analysis(B):
    x = []
    y = []
    for n in range(2,B):
        start = 10**(n)
        end = 10**(n+1)
        x.append(n)
        L = random.randint(start,end)
        strt = time.perf_counter()
        ans =  fermat(L,10)
        end = time.perf_counter()
        y.append(end-strt)
    plt.figure(figsize=(15,10))
    plt.plot(x,y, color = 'b')
    plt.title("Time variation for fermat(n,k)")
    plt.xlabel('No. of Digits in N')
    plt.ylabel('Time (in seconds)')
    plt.show()


# In[52]:


analysis(500)

