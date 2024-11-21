#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


C_to_K = 273.15
c_p_dry = 1005.7 # J/kgK
c_V_dry = 718.66 # J/kgK
eps = 0.6220
k_dry = 0.2854 


# In[11]:


def sat_vapor_pressure(T):
    '''
    T in C
    '''
    e_s = 6.112 * np.exp((17.67 * T)/(T + 243.5))
    return e_s # mb

def sat_vapor_temperature(e_s):
    '''
    e_s in mb
    '''
    T = ((243.5 * np.log(e_s)) - 440.8)/(19.48 - np.log(e_s))
    return T # C

def sat_mixing_ratio(p, T):
    '''
    p in mb
    T in C
    '''
    e_s = sat_vapor_pressure(T)
    w_s = eps * (e_s/(p - e_s))
    return w_s # in g/kg

def mixing_ratio_line(p, w_s):
    '''
    p in mb
    w_s is unitless (kg/kg)
    '''
    #w_s = w_s*1000
    #e_s = ((w_s * p)/eps)/(1 + (w_s/eps))
    e_s = (p * w_s)/(eps + w_s)
    T = sat_vapor_temperature(e_s)
    return T # C

def RH(T, p, w):
    '''
    T in C
    p in mb
    w_s is unitless (kg/kg)
    '''
    e = (p * w)/(eps + w)
    e_s = sat_vapor_pressure(T)
    RH = 100 * (e/e_s)
    return RH # %

def T_LCL(T, RH):
    '''
    T in K
    RH in %
    '''
    T_LCL = 55 + (1/((1/(T-55)) - (np.log(RH/100)/2840)))
    return T_LCL # K

def theta(T, p):
    '''
    T in C
    p in mb
    '''
    w_s = sat_mixing_ratio(p, T)
    Tk = T + C_to_K
    a = 0.2854 * (1 - (0.28 * (10**-3) * w_s))
    theta = Tk * np.power((1000/p), a)
    return theta
    
def theta_dry(theta, p, p_0=1000.0):
    theta_dry = theta * np.power(p_0/p, -k_dry)
    return theta_dry # K

def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    '''
    T in C
    p in mb
    w in kg/kg
    eq 43 in Bolton
    '''
    Tk = T + C_to_K
    RH1 = RH(T, p, w)
    T_LCL1 = T_LCL(Tk, RH1)
    #a = 0.2854 * (1 - (0.28 * (10**-3) * w))
    #theta = Tk * (p_0/p) ** a
    #theta_ep = theta * np.exp((3.376/T_LCL1)-0.00254) * w * (1 + (0.81 * (10**-3) * w))
    A = Tk*(p_0/p)**(0.2854*(1-(0.28*10**-3)*w))
    B = np.exp(((3.376/T_LCL1)-0.00254) * w*(1+(0.81*10**-3)*w))
    return A*B
    #return theta_ep # K


def theta_ep_field(T, p, p_0=1000.0):
    w_s = sat_mixing_ratio(p, T) * 1000
    print(w_s)
    theta_ep_field = pseudoeq_potential_T(T, p, w_s, p_0)
    print(theta_ep_field)
    return theta_ep_field

# In[ ]:




