#!/usr/bin/env python
# coding: utf-8

# ## Miller-Rabin test for compositeness

# In[55]:


import random 
def rabin_miller(N,k): # N : Any odd integer
        m=(N-1)/2
        h=1
        while m%2==0:
            m=int(m/2)
            h = h+1
        m = int(m)
        for i in range(k):
            x = random.randint(2,N-1)
            if pow(x,m,N)!=1:
                l=0
                for j in range(h):
                    if pow(x,pow(2,j)*m,N)!=N-1:
                        l=l+1
                if l==h:
                    return "Composite"
        return "Undecided"


# In[56]:


rabin_miller(100123456789,5)


# In[57]:


rabin_miller(232354763543765732524573251,15)


# In[59]:


import matplotlib.pyplot as plt
import time
import math
def analysis(B):
    x = []
    y = []
    for n in range(2,B):
        start = 10**(n)
        end = 10**(n+1)
        x.append(n)
        L = random.randint(start,end)
        while L%2 == 0:
            L = random.randint(start,end)
        strt = time.perf_counter()
        res = rabin_miller(L,5)
        end = time.perf_counter()
        y.append(end-strt)
    plt.figure(figsize=(15,10))
    plt.title("rabin_miller(N,k) Execution Time Variation for odd N")
    plt.plot(x,y, color = 'b')
    plt.xlabel('No. of Digits in N')
    plt.ylabel('Time Taken (in seconds)')
    plt.show()


# In[60]:


analysis(200)

