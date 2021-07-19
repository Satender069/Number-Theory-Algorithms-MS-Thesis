#!/usr/bin/env python
# coding: utf-8

# # **Quadratic Sieve Algorithm**
# 

# For number $n$ to be factored it aims to find the solution to $x^2\equiv y^2\ (mod\ n)$ $\implies$ $(x-y)(x+y)=0\ (mod\ n)$. If so, we can compute $gcd(x-y,n)$ and $gcd(x+y,n)$ to see if these are non-trivial divisors. 
# 
# First, we define the quadratic polynomial $A(x) = x^2-n$ to work with and compute $A(x_1),A(x_2),....,A(x_k)$ for some $x_i's$, that we will define below, and pick up a subset $s.t.$ $A(x_{i_1})A(x_{i_2})....A(x_{i_r})$ is a square, say, $y^2$. From our equation of $A(x)$ we then get: 
# 
# $$
# y^2=A(x_{i_1})A(x_{i_2})....A(x_{i_r}) \equiv (x_{i_1}x_{i_2}....x_{i_r})^2 
# $$
# 
# Hence, we check for $gcd(y-(x_{i_1}x_{i_2}....x_{i_r}),n)$ to be non-trivial. If not, then we look for different subset of $A(x_i)$ and repeat the process.

# ### Setting up a Sieving Interval and a Factor Base
# 
# We need efficient way to determine our $x_i$ as to get a perfect square as we are going to factor $f(x_i)$ with the primes in our factor base so we want it to be small and to factor over set of small prime numbers. 
# 
# To make $f(x_i)= {x_i}^2-n$ small we use $x_i = ( \lfloor \sqrt n \rfloor +r )$ where $r\in [-M,M]$, which is our $sieving\ interval$. 
# 
# For negative $f(x_i)$ values we have to include $-1$ in our factor base. Here we take care of sign by calculating absolute value of $f(x_i)'s$.
# 
# Now, for $factor\ base$ we need primes $p$, below some bound $B$, $s.t.$ $n$ is a quadratic residue modulo $p$ as others do not divide any of the $f(x_i)'s$. Hence, we need legendre symbol $\big(\frac{n}{p}\big)$=1.
# 
# 
# 
# 

# ### Finding suitable list of primes upto bound B

# In[1]:


def jacobi(p,a):   # to check whether 'a' is a  Quadratic Residue modulo p 
        b=a%p
        c=p
        s=1
        while b>=2:
            while b%4==0:
                b = b/4
            if b%2 == 0:
                if c%8==3 or c%8==5:
                    s = -s
                b=b/2
            if b==1:
                break
            if b%4==c%4==3:
                s=-s
            b1=b
            b=c%b
            c=b1
        return s


# Below program is to find primes upto $B$ by trial method. While working with large factor base one should replace it with better versions for the same. 

# In[2]:


def primes(B,n):       
    factor_base = [2]
    for num in range(3,B + 1):
        for i in range(2, num):
            if (num % i) == 0:
                break
        else:
            if jacobi(num,n) == 1:
                factor_base.append(num)
    return factor_base


# ### Working with Sieving interval [-M,M]

# After getting out $factor\ base$ we now need to calculate $f(x_i)$ values for $( \lfloor \sqrt n \rfloor - M )\leq x_i\leq ( \lfloor \sqrt n \rfloor + M)$. We keep only those values which are factored within the primes of our factor base only.
# 
# Assume our factor base to be $(a_1,a_2,.....,a_k)$ and factors factoring within these primes be:
# 
# $$
# f(x_1) = (a_1)^{r_{11}} (a_2)^{r_{12}}.... (a_k)^{r_{1k}} \\
# f(x_2) = (a_1)^{r_{21}} (a_2)^{r_{22}}.... (a_k)^{r_{2k}} \\
# . \\
# . \\
# . \\
# f(x_p) = (a_1)^{r_{p1}} (a_2)^{r_{p2}}.... (a_k)^{r_{pk}}
# $$
# 
# Below we return lists of the powers corresponding to prime factors for each $f(x_i)$ and,  all $x_i$  and their $f(x_i)$ values.

# In[3]:


import math
def factoring(B,M,n): # B- Bound for prime factor base, M- sieving interval bound, n- number to be factored
    power = 0
    dic = {}
    list_final=[]
    fx_list=[]
    x_list=[]
    gif = math.floor(pow(n,1/2))
    factor_base = primes(B,n)
    for i in range(gif-M,gif+M+1):
        num = i*i-n
        for j in factor_base:            # collecting powers of primes of factor base
            while num%j == 0:
                power += 1
                num = num/j
            dic["{}".format(j)]= power
            power = 0
            num = i*i-n
        
        number = 1
        for j in factor_base:             # checking if number factored within our factor base or not
            number = number*pow(j,dic["{}".format(j)])
        if number == abs(num):
            x_list.append(i)
            fx_list.append(abs(num))
            list_final.append(dic.copy())
    return list_final,x_list,fx_list


# In[4]:


factoring(30,16,9487)


# In[5]:


from sympy import primefactors as pf
pf(8787)


# In[6]:


import math as m
m.gcd(94+7,8787)


# After having factored required $f(x_i)'s$ we need to check for the subset, say $S = \{f(x_{i_1}),f(x_{i_2}),....f(x_{i_s})\}$ of these values $s.t.$  $s\ =\ f(x_{i_1})f(x_{i_2})....f(x_{i_s})$ is a perfect square. Which means we need the powers of each prime in $s$ to be even. 
# 
# So, let $M$ be our matrix of powers of all factored $f(x_i)'s$ as below:
# 
# $\begin{bmatrix}
# r_{11}& r_{12} & r_{13} & ..... & r_{1k}\\
# r_{21}& r_{22} & r_{23} & ..... & r_{2k}\\
# . & . & . & & .\\
# . & . & . & & .\\
# r_{p1}& r_{p2} & r_{p3} & ..... & r_{pk}\\
# \end{bmatrix}$
# 
# And let $M'\ =\ M\ (mod\ 2)$ be the matrix with entries $r_{ij}\ (mod\ 2);\ 1\leq i\leq p,\ 1\leq j\leq k$. 

# ##### Gaussian Elimination
# 
# We will now use gaussian elimination to get desired subset $S$ out of $f(x_i)'s$. For this we will work with matrix $M'$ and identity matrix $A_{p\times p}$ modulo 2. Each row vector of $A$ consists of $0's$ and $1's$. We use $A$ to keep track of the combination of $f(x_i)'s$ that gives us a perfect square.
# 
# Starting with first column, we find the first row with $1$ and then add modulo $2$ to succeeding rows with $1's$ in that column. We do the same with corresponding rows of matrix $A$ which keeps the track of which $f(x_i)'s$ are multiplied. Then we eliminate that row from the matrix $M'$.
# 
# Then we repeat the above step for column $2$. Notice that after going through each column the remaining rows have only $0$ entries in that column.
# 
# We keep repeating the above process until we get a row with all the zero entries, which would eventually happen if we take more $f(x_i)$ values than in factor base. Hence the bounds are very important for quadratic sieve algorithm. Zero entries imply that the power of primes is even.
# 
# Assume $j^{th}$ to be the required row then we check $j^{th}$ row's entries in $A$ for $f(x_i)$ values which on multiplying will give us the required square. Now we can solve for $x^2\equiv y^2\ (mod\ n)$.
# 
# 
# If we do not get a non-trivial factor then we repeat the process further we find another such row.
# 
# 
# 
# 

# In[7]:


import numpy as np
def qsa(B,M,n):                                  # Bound for factor base, Bound for sieving interval, odd integer to be factored
    list1 = []
    final_list,x_values,fx_values = factoring(B,M,n) # vectors list, x-values, f(x) values
    factor_base = primes(B,n)
    
    len_final_list = len(final_list)                 # no. of rows
    len_factor_base= len(factor_base)                # factor_base length
    
    #....... Setting Matrices ..............
    
    for i in range(len_final_list):
        list1.append(list(final_list[i].values()))
    main_matrix = np.array(list1)                    # Matrix with vectors corresponding to factors powers before modulo 2
    #................................................
    working_matrix = (main_matrix.copy()) % 2        # Matrix After modulo 2 calculation
    tracking_matrix = np.identity(len_final_list)    # Identity matrix for row operations history
    
    
    worked_row_list = []                             # List of row numbers that have been eliminated.
    for column_no in range(len_factor_base):
        for row_no in range(len_final_list):
            if row_no not in worked_row_list:
                
                curr_row = working_matrix[row_no,]
                
                for row in range(row_no,len_final_list):
                    if row not in worked_row_list:
                        if all([i==0 for i in (working_matrix[row,]%2)]): 
                            
                            # .............FINAL X VALUE..............
                            x_value = 1
                            for i,j in zip(tracking_matrix[row,],x_values):
                                if i!= 0:
                                    x_value = (x_value*j)%n
                                    
                            #.............FINAL Y VALUE............
                            y_value = 1
                            sum_rows=np.zeros((1,len_factor_base))
                            for i in range(len_final_list):
                                if tracking_matrix[row,i]==1:
                                    sum_rows += main_matrix[i,]
                            for power,factor in zip(sum_rows[0],factor_base):
                                y_value *= pow(factor,int(power/2),n)
                                     
                            # .... Checking if Non-trivial Divisor......       
                            divisor =  math.gcd(x_value+y_value,n) 
                            if divisor not in [1,n]:
                                return "Divisors of {} are = {} and {} ".format(n,divisor,int(n/divisor))
                                break
                            
                #....... Gaussian Elimination..............
                if working_matrix[row_no,column_no] == 1:
                    for r in range(row_no+1,len_final_list):
                        if working_matrix[r,column_no]==1:
                            working_matrix[r,]=(working_matrix[r,]+curr_row)%2
                            tracking_matrix[r,]=(tracking_matrix[r,]+ tracking_matrix[row_no,])%2
                    worked_row_list.append(row_no)


# In[8]:


import time
strt = time.perf_counter()
print(qsa(10,10,4295229443))
end = time.perf_counter()
print(f'Time taken = {end-strt}')


# In[9]:


import time
strt = time.perf_counter()
print(qsa(80,1720,1234567895341))
end = time.perf_counter()
print(f'Time taken = {end-strt}')


# In[40]:


print(qsa(1000,100000,1000000000099987889))

