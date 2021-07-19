#!/usr/bin/env python
# coding: utf-8

# ### Example : Find  irreducible factors of $f(x) = 8x^{10} + 28x^9 + 24x^8 + 36x^7 + 14x^6 + 35x^5 - 12x^4 + 22x^3 - 13x^2 + 5x - 3$ in $\mathbb{Z}[X]$.
# $\textit{Solution:}$ Factors are: $(x + 3),(2x - 1),(2x^2 + 1),(x^4 + x^3 + x^2 + 1)$
# 

# In[34]:


f = Poly(8*x**10 + 28*x**9 + 24*x**8 + 36*x**7 + 14*x**6 + 35*x**5 - 12*x**4 + 22*x**3 - 13*x**2 + 5*x - 3, x, domain='ZZ')


# ##### Note: $f$ is a global variable.

# In[35]:


import sympy as sp
from sympy import prime 
from sympy.abc import x
from sympy import *
from sympy import GF
from sympy import sqf_part
from sympy import factor_list


# ##### Getting square-free part of $f(x)$.

# In[36]:


f = sqf_part(f)
f


# #### 0). Choosing suitable $'p'$ values such that $f(x)\ (mod p)$ is square-free in $\mathbb{Z}_p[X]$.

# In[52]:


def factor_base(f, pbound):
    factor_base = []
    poly_list = []
    prime_list = [prime(i) for i in range(1,pbound+1)]
    for num in prime_list:
        coeff_f = f.all_coeffs()
        
        if coeff_f[0]%num != 0:
            
            #Reducing f(x) mod(num)
            K = GF(num)                         
            g = sp.Poly(coeff_f,x,domain = K)
        
            # Checking if square_free
            if sp.discriminant(g)!= 0:
                poly_list.append(g)
                factor_base.append(num)
    return factor_base,poly_list      


# In[38]:


factor_base(f,10)


# ### Factoring $f(x)\ mod(p)$ in $\mathbb{Z}_p[X]$.

# #### Distinct-degree Factorization(DDF)

# In[39]:


from sympy import monic
from sympy import ZZ


def DDF(f,p):  # Enter 'square-free', 'primitive' f(x) \in Z[X].
    
    K = GF(p) 
    f_modp = Poly(f, domain = K)
    f = monic(f_modp,domain=K)   # Converting to Monic in K.
    f_coef = f_modp.all_coeffs()
    degree_f = f.degree()
    
    # We have monic,square-free POLYNOMIAL =  f, COEFFICIENTS = f_coef, PRIME = p
    F = f
    DDF_list = []
    degree = 1
    deg_list = []
    total_deg = 0
    while total_deg != degree_f:
        
        if p**degree > f.degree():
            B = sp.poly(x**p)
            for i in range(1,degree):
                B = sp.rem(B**p,f,domain = K)
        else:
            B = sp.poly(x**(p**degree),domain = K)
        factor = sp.gcd(f,Poly(B-sp.Poly([1,0],x)),domain = K)
        if factor != 1:
            deg_list.append(degree)
            DDF_list.append(factor)
            total_deg = total_deg + degree
        degree = degree+1
        
        f = sp.quo(f,factor)
        if f == 1:
            break
    return DDF_list, deg_list,p


# #### Equal-Degree Factorization(EDF)

# In[40]:


import numpy as np
from numpy.polynomial import polynomial as P
from numpy import random as r


def EDF(f,p):
    
    A_d,deg_list,p = DDF(f,p)
    factors_list = []
    irr_factors_list = []
    K = GF(p)
    n_factors = 0
    if p>2: 
        
        for (poly,deg) in zip(A_d,deg_list):
            
            # Checking if Poly is an irreducible polynomial 
            if deg == poly.degree():
                irr_factors_list.append(poly)
                n_factors = n_factors +1
            else:
                n_factors = n_factors + (poly.degree()/deg)
                T = sp.random_poly(x,poly.degree(),inf = 0,sup = p,domain = K)
                d = 2
                d1 = ((p**deg)-1)/2
                T_2 = T**2
                while d <= d1:
                    T_2 = sp.rem(T_2,poly)*T
                    d = d+1
                F1 = sp.gcd(T_2+1,poly,domain = K)
                F2 = sp.gcd(T_2-1,poly,domain = K)
                F3 = sp.quo(poly,F1*F2,domain = K)
                
                deg_f1 = F1.degree()
                deg_f2 = F2.degree()
                deg_f3 = F3.degree()
                
                if deg_f1 != 0:
                    if deg_f1 == deg:
                        irr_factors_list.append(F1)
                    else:
                        factors_list.append([F1,deg])
                if deg_f2 != 0:
                    if deg_f2 == deg:
                        irr_factors_list.append(F2)
                    else:
                        factors_list.append([F2,deg])
                if deg_f3 != 0:
                    if deg_f3 == deg:
                        irr_factors_list.append(F3)
                    else:
                        factors_list.append([F3,deg])
    return factors_list,irr_factors_list,p, n_factors


# #### Factoring completely using DDF and EDF as above

# ## ======================== =========== ==================================

# Programs below are to be hadled for each prime $p$ separately.

# In[53]:


factors_list, irr_factors_list,p,n_factors = EDF(f,5)  # Input : polynomial in Z[X], primes bound, prime number
K = GF(p)
while len(irr_factors_list) != n_factors:
    list2 = []
    for elt in factors_list:
        polynomial = elt[0]
        degree = elt[1]
        T = sp.random_poly(x,degree,inf = 0,sup = p,domain = K)
        d = 2
        d1 = ((p**degree)-1)/2
        T_2 = T**2
        while d <= d1:
            T_2 = sp.rem(T_2,polynomial)*T
            d = d+1
        F1 = sp.gcd(T_2+1,polynomial,domain = K)
        F2 = sp.gcd(T_2-1,polynomial,domain = K)
        F3 = sp.quo(polynomial,F1*F2,domain = K)       
        deg_f1 = F1.degree()
        deg_f2 = F2.degree()
        deg_f3 = F3.degree()
        if deg_f1 != 0:
            if deg_f1 == degree:
                irr_factors_list.append(F1)
            else:
                list2.append([F1,degree])
        if deg_f2 != 0:
            if deg_f2 == degree:
                irr_factors_list.append(F2)
            else:
                list2.append([F2,degree])
        if deg_f3 != 0:
            if deg_f3 == degree:
                irr_factors_list.append(F3)
            else:
                list2.append([F3,degree])
    factors_list = list2.copy()
print("FACTORS=", irr_factors_list)


# ##### Storing all possible (h modp) in list H

# In[42]:


H = []
irr_factors_list
for factor in irr_factors_list:
    if sp.rem(f,Poly(factor,domain = ZZ)) != 0:
        H.append(factor)
H


# ##### Irreducible factors so far, (in Z[X]). Stored in irr_factors_list.

# In[43]:


for poly in H:
    irr_factors_list.remove(poly)
irr_factors_list


# ###### Note: Lists H and irr_factors_list are global variables.

# # =============================================================

# #### Calculating k for suitable p^k value

# In[54]:


import scipy.special as ss
from sympy import Matrix
from sympy import ceiling,Float
def start(f,H,r,p):  # Input : f - polynomial, H[r]- Entries from H, p - prime number from factor_base.
    h = H[r]
    deg_h = h.degree()
    n = f.degree()
    m = n-1
    M = Matrix(list(f.coeffs()))
    B = 2**(m*n/2) * (ss.binom(2*m,m))**(n/2) * M.norm()**(m+n)
    k = ceiling(ceiling(sp.log(B+1,p))/deg_h) 
    return k


# #### Hensel lifting 

# #### Lifting $(mod p)$ to $(mod p^2)$.

# In[55]:


def hensel_1(f,g,p): # Input: f - Polynomial, g - Monic factor in Z_p[X], p - Prime from factor base
    K = GF(p)
    K_1 = GF(p*p)
    h = Poly(sp.quo(f,g,domain = K),domain = ZZ)
    t,s,l = sp.gcdex(g,h,domain = K)
    g = Poly(g,domain = ZZ)
    e = Poly(f-g*h,domain = K_1)
    q,r = sp.div(s*e,g,domain = K_1)
    g1 = Poly(g+r,domain = K_1)
    h1 = Poly(h+t*e+q*h,domain = K_1)
    return g1,h1


# #### Lifting $(mod q)$ to $(mod qr)$ where $r = gcd(p,q)$ for some $p,q$ ( not necessarily prime).

# In[56]:


import math as m
def hensel_2(u,v,q,p):  # u = vw in Z_q s.t. av+bw = 1 in Z_p
    r = m.gcd(q,p)
    K = GF(q)
    K_1 = GF(q*r)
    w = Poly(sp.quo(u,v,domain = K),domain = ZZ)
    a,b,l = sp.gcdex(v,w,domain = GF(p))
    b = Poly(b,domain = ZZ)
    v = Poly(v,domain = ZZ)
    f = Poly(sp.quo(u-v*w,sp.Poly([q],x),domain = ZZ))
    t= sp.quo(b*f,v) 
    v_1 = b*f - t*v
    w_1 = a*f + t*w
    V = v + q*v_1
    W = w + q*w_1
    V = Poly(V,domain = K_1)
    W = Poly(W,domain = K_1)
    return V,W


# In[47]:


def hmodpk(f,H,r,p):  # Input : f,H[r],p as in above programs

    k = start(f,H,r,p)
    h = Poly(H[r],domain = ZZ)
    o = 1
    while 2**o < k:  # k factored as p^(2^(r-1)).p^(k-2^(r-1))
        o = o+1
    n = 1
    q = p
    while n != o:
        h_1,g_1 = hensel_1(f,h,q)
        q = q*q
        h = Poly(h_1,domain = ZZ)
        n = n+1
    h_1,g_1 = hensel_2(f,h,q,int(p**(k-2**(o-1)))   )
    return h_1


# ###  Getting irreducible factor in $\mathbb{Z}[X]$ for each suitable $h\ (mod\ p)$

# In[48]:


import olll

def red_lattice(f,H,r,p): 
    
    h_mod_pk = hmodpk(f,H,r,p)
    n = f.degree()
    l = h_mod_pk.degree()
    
    # Creating set of m values
    u = 0
    M = []
    while l <= (n-1)/(2**u):
        M.append(int((n-1)/(2**u)))
        u = u+1
    M.reverse()
    k = start(f,H,r,p)
    f_norm = Matrix(list(f.all_coeffs())).norm()
    
    for m in M:
        
        # CREATING LATTICE BASIS VECTORS
        P1 = []
        P2 = []
        for i in range(l):
            if i != 0:
                P1.append((p**k)*sp.poly(x**i))
            else:
                P1.append((p**k)*sp.Poly([1],x))
        for j in range(m-l+1):
            if j != 0:
                P2.append(Poly(h_mod_pk,domain = ZZ)*sp.poly(x**j))
            else:
                P2.append(Poly(h_mod_pk,domain = ZZ)*sp.Poly([1],x)) 
                
        # Managing dimension of vectors and collecting in one list                
        Lattice = []             
        dim_vectors = m+1 
        for poly in P1:
            poly_list = poly.all_coeffs()
            poly_list.reverse()
            for i in range(m+1):
                if i > len(poly_list)-1:
                    poly_list.append(0)
            Lattice.append(poly_list)
        for poly in P2:
            poly_list = poly.all_coeffs()
            poly_list.reverse()
            for i in range(m+1):
                if i > len(poly_list)-1:
                    poly_list.append(0)
            Lattice.append(poly_list)
            
        # GETTING REDUCED LATTICE BASIS VECTORS
        Lattice2 = []
        for vec in Lattice:
            L = []
            for num in vec:
                L.append(int(num))
            Lattice2.append(L)
        reduced_lattice = olll.reduction(Lattice2,3/4)
        
        # CHECKING IF deg(h0) <= m
        lattice_poly = []
        bound = (p**(k*l)/(f_norm**m))**(1/n)

        t = 0
        for vector in reduced_lattice:

            if Matrix(vector).norm() < bound:
                vector = list(vector)
                vector.reverse()
                lattice_poly.append(sp.Poly(vector,x,domain = ZZ))
                t = t+1
            else:
                break
        if t != 0:
            H_0 = lattice_poly[0]
            for i in range(t):
                H_0 = sp.gcd(H_0,lattice_poly[i])
            return H_0


# ### Factoring completely into irreducibles in $\mathbb{Z}[X]$

# In[49]:


def factorization(f,H,p):
    h_0 = red_lattice(f,H,0,p)
    print(h_0)
    fact_check = [h_0]
    for i in range(1,len(H)):
        if sp.rem(h_0,H[i]) != 0:
            h_0 = red_lattice(f,H,i,p)
            if h_0 not in fact_check:
                print(h_0)


# In[50]:


factorization(f,H,5)
print(irr_factors_list)

