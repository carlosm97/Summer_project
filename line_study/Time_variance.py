#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 17:51:16 2021

@author: charlie
"""

import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
import math as mth
import os
os.chdir('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/pyIACOB')
from scipy.optimize import curve_fit
from scipy import integrate
from db import *
from spec import *
from RV import *
from line_study import *
import progressbar
import time

# In[Telescope_effect_SiIII]


SM=reading_summary('../stars/HD2905/4574/00_HD2905_4574.txt')
m_SiIII_fm,f_SiIII_fm,s_SiIII_fm=[],[],[]
m_HJD,f_HJD,s_HJD=[],[],[]
error_m,error_f,error_s=[],[],[]
for i in range(len(SM.first_moment)):
    if SM.Instrument[i]=="'MERCATOR'":
        m_SiIII_fm.append(SM.first_moment[i]),m_HJD.append(SM.HJD[i]),error_m.append(SM.error_fm[i])
    elif SM.Instrument[i]=="'FIES'":
        f_SiIII_fm.append(SM.first_moment[i]),f_HJD.append(SM.HJD[i]),error_f.append(SM.error_fm[i])
    elif SM.Instrument[i]=="'SONG'":
        s_SiIII_fm.append(SM.first_moment[i]),s_HJD.append(SM.HJD[i]),error_s.append(SM.error_fm[i])
    else:
        print('Unknown instrument')
plt.figure()
#plt.plot(np.vstack(m_HJD)-2.45e6,m_SiIII_fm,'xr',label='HERMES')        
#plt.plot(np.vstack(f_HJD)-2.45e6,f_SiIII_fm,'xb',label='FIES')            
#plt.plot(np.vstack(s_HJD)-2.45e6,s_SiIII_fm,'xg',label='SONG')

plt.plot(np.vstack(m_HJD)-2.45e6,np.vstack(m_SiIII_fm)-np.mean(m_SiIII_fm),'.r',label='Rescaled HERMES')        
plt.plot(np.vstack(f_HJD)-2.45e6,np.vstack(f_SiIII_fm)-np.mean(f_SiIII_fm),'.b',label='Rescaled FIES')            
plt.plot(np.vstack(s_HJD)-2.45e6,np.vstack(s_SiIII_fm)-np.mean(s_SiIII_fm),'.g',label='Rescaled SONG')
plt.grid(),plt.legend()

plt.ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
plt.xlabel('HJD [days]-2450000')

Z_m=sorted(zip(np.vstack(m_HJD)-2.45e6,np.vstack(m_SiIII_fm)-np.mean(m_SiIII_fm)))
Z_f=sorted(zip(np.vstack(f_HJD)-2.45e6,np.vstack(f_SiIII_fm)-np.mean(f_SiIII_fm)))
Z_s=sorted(zip(np.vstack(s_HJD)-2.45e6,np.vstack(s_SiIII_fm)-np.mean(s_SiIII_fm)))

m_HJD,m_fm=[Z_m[i][0] for i in range(len(Z_m))],[Z_m[i][1] for i in range(len(Z_m))]
f_HJD,f_fm=[Z_f[i][0] for i in range(len(Z_f))],[Z_f[i][1] for i in range(len(Z_f))]
s_HJD,s_fm=[Z_s[i][0] for i in range(len(Z_s))],[Z_s[i][1] for i in range(len(Z_s))]

# Let's use this Si III line to study the variance in different ways. 

# In[random]
#plt.close('all')
# Consider the std with changing random data. 
np.random.seed(5)
STD_s_rn,STD_f_rn,STD_m_rn=[],[],[]
fig, axs=plt.subplots(sharex=True)
for i in range(50000):
    # Selection of the data to be used:
    num=[]
    rnd=np.random.randint(0,len(s_HJD),int(0.1*len(s_HJD)))
    for j in rnd:
        num.append(s_fm[j])
    STD_s_rn.append(np.std(num))
    num=[]
    rnd=np.random.randint(0,len(f_HJD),int(0.1*len(f_HJD)))
    for j in rnd:
        num.append(f_fm[j])
    STD_f_rn.append(np.std(num))
    num=[]
    rnd=np.random.randint(0,len(m_HJD),int(0.2*len(m_HJD)))   
    for j in rnd:
        num.append(m_fm[j])
    STD_m_rn.append(np.std(num))    
    
axs.hist(STD_s_rn,bins=50,density=True,alpha=0.3,label='SONG')

axs.hist(STD_f_rn,bins=50,density=True,alpha=0.3,label='FIES')

axs.hist(STD_m_rn,bins=50,density=True,alpha=0.3,label='HERMES')
axs.legend()
# In[convolution]


STD_s_cv,STD_f_cv,STD_m_cv=[],[],[]
#step = 40
data_0=min(f_HJD[0],m_HJD[0],s_HJD[0])[0]
data_1=max(f_HJD[-1],m_HJD[-1],s_HJD[-1])[0]
#n_step=int((data_1-data_0)/step)+1
f_HJD,m_HJD,s_HJD,f_fm,m_fm,s_fm=np.vstack(f_HJD),np.vstack(m_HJD),np.vstack(s_HJD),np.vstack(f_fm),np.vstack(m_fm),np.vstack(s_fm)

for n in range(1,50):
    STD_f_cv.append(np.nanstd(f_fm[int(0.02*n*len(f_fm)):int((0.02*n+0.1)*n*len(f_fm))]))
    STD_m_cv.append(np.nanstd(m_fm[int(0.02*n*len(m_fm)):int((0.02*n+0.1)*len(m_fm))]))
    STD_s_cv.append(np.nanstd(s_fm[int(0.02*n*len(s_fm)):int((0.02*n+0.1)*n*len(s_fm))])) 
    # Buscaremos que sea mayor que data_0+3*n y menor que data_0+3*(n+1), sumaremos y haremos nanstd, e histograma 
'''
try: plt.close(hist2)
except: pass
hist2=plt.figure()
plt.hist(STD_m_cv,density=True,alpha=0.3,label='HERMES') 
plt.hist(STD_f_cv,density=True,alpha=0.3,label='FIES') 
plt.hist(STD_s_cv,density=True,alpha=0.3,label='SONG') 
plt.legend()
'''
try: plt.close(sf)
except: pass
sf,sa=plt.subplots(2,sharex=True)

sa[0].hist(STD_s_rn,bins=50,density=True,alpha=0.3,label='SONG')
sa[0].hist(STD_f_rn,bins=50,density=True,alpha=0.3,label='FIES')
sa[0].hist(STD_m_rn,bins=50,density=True,alpha=0.3,label='HERMES')
sa[0].legend()

sa[1].hist(STD_s_cv,bins=20,density=True,alpha=0.3,label='SONG')
sa[1].hist(STD_f_cv,bins=20,density=True,alpha=0.3,label='FIES') 
sa[1].hist(STD_m_cv,bins=20,density=True,alpha=0.3,label='HERMES') 

sa[1].legend()
plt.show()
# In[convolution_2_sizes]
STD_s_1,STD_f_1,STD_m_1=[],[],[]
STD_s_10,STD_f_10,STD_m_10=[],[],[]
STD_s_25,STD_f_25,STD_m_25=[],[],[]
step1,step10,step25 = 1,10,25
data_0=min(f_HJD[0],m_HJD[0],s_HJD[0])
data_1=max(f_HJD[-1],m_HJD[-1],s_HJD[-1])
n_step_1,n_step_10,n_step_25=int((data_1-data_0)/step1)+1,int((data_1-data_0)/step10)+1,int((data_1-data_0)/step25)+1
for i in range(n_step_1):
    mask1_s=np.logical_and(np.vstack(s_HJD)>(i+s_HJD[0]),np.vstack(s_HJD)<(i+step1+s_HJD[0]))
    if len(np.vstack(s_fm)[mask1_s])>1:
        STD_s_1.append(np.std(np.vstack(s_fm)[mask1_s]))
    mask1_f=np.logical_and(np.vstack(f_HJD)>(i+f_HJD[0]),np.vstack(f_HJD)<(i+step1+f_HJD[0]))
    if len(np.vstack(f_fm)[mask1_f])>1:
        STD_f_1.append(np.std(np.vstack(f_fm)[mask1_f]))
    mask1_m=np.logical_and(np.vstack(m_HJD)>(i+m_HJD[0]),np.vstack(m_HJD)<(i+step1+m_HJD[0]))
    if len(np.vstack(m_fm)[mask1_m])>1:
        STD_m_1.append(np.std(np.vstack(m_fm)[mask1_m]))

for i in range(n_step_10):
    mask10_s=np.logical_and(np.vstack(s_HJD)>(step10*i+s_HJD[0]),np.vstack(s_HJD)<(step10*(i+1)+s_HJD[0]))
    if len(np.vstack(s_fm)[mask10_s])>1:
        STD_s_10.append(np.std(np.vstack(s_fm)[mask10_s]))
    mask10_f=np.logical_and(np.vstack(f_HJD)>(step10*i+f_HJD[0]),np.vstack(f_HJD)<(step10*(i+1)+f_HJD[0]))
    if len(np.vstack(f_fm)[mask10_f])>1:
        STD_f_10.append(np.std(np.vstack(f_fm)[mask10_f]))
    mask10_m=np.logical_and(np.vstack(m_HJD)>(step10*i+m_HJD[0]),np.vstack(m_HJD)<(step10*(i+1)+m_HJD[0]))
    if len(np.vstack(m_fm)[mask10_m])>1:
        STD_m_10.append(np.std(np.vstack(m_fm)[mask10_m]))

for i in range(n_step_25):
    mask25_s=np.logical_and(np.vstack(s_HJD)>(step25*i+s_HJD[0]),np.vstack(s_HJD)<(step25*(i+1)+s_HJD[0]))
    if len(np.vstack(s_fm)[mask25_s])>1:
        STD_s_25.append(np.std(np.vstack(s_fm)[mask25_s]))
    mask25_f=np.logical_and(np.vstack(f_HJD)>(step25*i+f_HJD[0]),np.vstack(f_HJD)<(step25*(i+1)+f_HJD[0]))
    if len(np.vstack(f_fm)[mask25_f])>1:
        STD_f_25.append(np.std(np.vstack(f_fm)[mask25_f]))
    mask25_m=np.logical_and(np.vstack(m_HJD)>(step25*i+m_HJD[0]),np.vstack(m_HJD)<(step25*(i+1)+m_HJD[0]))
    if len(np.vstack(m_fm)[mask25_m])>1:
        STD_m_25.append(np.std(np.vstack(m_fm)[mask25_m]))
        
sp,sa=plt.subplots(4,sharex=True)
'''
sa[0].hist(STD_s_1,bins=20,density=True,alpha=0.3,label='SONG')
sa[0].hist(STD_f_1,bins=20,density=True,alpha=0.3,label='FIES')
sa[0].hist(STD_m_1,bins=20,density=True,alpha=0.3,label='HERMES')
sa[0].legend()
sa[1].set_title('10 days step')
sa[1].hist(STD_s_10,bins=20,density=True,alpha=0.3,label='SONG')
sa[1].hist(STD_f_10,bins=20,density=True,alpha=0.3,label='FIES')
sa[1].hist(STD_m_10,bins=20,density=True,alpha=0.3,label='HERMES')
sa[1].legend()
sa[2].set_title('25 days step')
sa[2].hist(STD_s_25,bins=20,density=True,alpha=0.3,label='SONG')
sa[2].hist(STD_f_25,bins=20,density=True,alpha=0.3,label='FIES')
sa[2].hist(STD_m_25,bins=20,density=True,alpha=0.3,label='HERMES')
sa[2].legend()
'''
maxi,mini=max(max(STD_s_1),max(STD_f_1),max(STD_m_1)),min(min(STD_s_1),min(STD_f_1),min(STD_m_1))
sa[0].set_title('1 day step')
H_s_1, bins_s_1 = np.histogram(STD_s_1,bins=20,range=[mini,maxi])
H_f_1, bins_f_1 = np.histogram(STD_f_1,bins=20,range=[mini,maxi])
H_m_1, bins_m_1 = np.histogram(STD_m_1,bins=20,range=[mini,maxi])
sa[0].bar((bins_s_1[1:]+bins_s_1[:-1])/2,H_s_1/max(H_s_1),alpha=0.3,width=np.mean(bins_s_1[1:]-bins_s_1[:-1]),label='SONG')
sa[0].bar((bins_f_1[1:]+bins_f_1[:-1])/2,H_f_1/max(H_f_1),alpha=0.3,width=np.mean(bins_f_1[1:]-bins_f_1[:-1]),label='FIES')
sa[0].bar((bins_m_1[1:]+bins_m_1[:-1])/2,H_m_1/max(H_m_1),alpha=0.3,width=np.mean(bins_m_1[1:]-bins_m_1[:-1]),label='HERMES')
sa[0].legend()


sa[1].set_title('10 day step')
H_s_10, bins_s_10 = np.histogram(STD_s_10,bins=20,range=[mini,maxi])
H_f_10, bins_f_10 = np.histogram(STD_f_10,bins=20,range=[mini,maxi])
H_m_10, bins_m_10 = np.histogram(STD_m_10,bins=20,range=[mini,maxi])
sa[1].bar((bins_s_10[1:]+bins_s_10[:-1])/2,H_s_10/max(H_s_10),alpha=0.3,width=np.mean(bins_s_10[1:]-bins_s_10[:-1]),label='SONG')
sa[1].bar((bins_s_10[1:]+bins_s_10[:-1])/2,H_f_10/max(H_f_10),alpha=0.3,width=np.mean(bins_f_10[1:]-bins_f_10[:-1]),label='FIES')
sa[1].bar((bins_s_10[1:]+bins_s_10[:-1])/2,H_m_10/max(H_m_10),alpha=0.3,width=np.mean(bins_m_10[1:]-bins_m_10[:-1]),label='HERMES')
sa[1].legend()


sa[2].set_title('25 day step')
H_s_25, bins_s_25 = np.histogram(STD_s_25,bins=20,range=[mini,maxi])
H_f_25, bins_f_25 = np.histogram(STD_f_25,bins=20,range=[mini,maxi])
H_m_25, bins_m_25 = np.histogram(STD_m_25,bins=20,range=[mini,maxi])
sa[2].bar((bins_s_25[1:]+bins_s_25[:-1])/2,H_s_25/max(H_s_25),alpha=0.3,width=np.mean(bins_s_25[1:]-bins_s_25[:-1]),label='SONG')
sa[2].bar((bins_s_25[1:]+bins_s_25[:-1])/2,H_f_25/max(H_f_25),alpha=0.3,width=np.mean(bins_f_25[1:]-bins_f_25[:-1]),label='FIES')
sa[2].bar((bins_s_25[1:]+bins_s_25[:-1])/2,H_m_25/max(H_m_25),alpha=0.3,width=np.mean(bins_m_25[1:]-bins_m_25[:-1]),label='HERMES')
sa[2].legend()



sa[3].set_title('Random numbers')
sa[3].hist(STD_s_rn,bins=50,density=True,alpha=0.3,label='SONG')
sa[3].hist(STD_f_rn,bins=50,density=True,alpha=0.3,label='FIES')
sa[3].hist(STD_m_rn,bins=50,density=True,alpha=0.3,label='HERMES')
sa[3].legend()

sa[3].set_xlabel(r'$\sigma_v$ [km/s]')
# In[prueba]

H, bins = np.histogram(STD_s_1,bins=20)
plt.figure()



# In[peak_peak]

np.random.seed(5)
pp_f,pp_m,pp_s=[],[],[]
for i in range(25000):
    # Selection of the data to be used:
    num=[]
    rnd=np.random.randint(0,len(s_HJD),int(0.1*len(s_HJD)))
    for j in rnd:
        num.append(float(s_fm[j]))
    pp_s.append(max(num)-min(num))
    num=[]
    rnd=np.random.randint(0,len(f_HJD),int(0.1*len(f_HJD)))
    for j in rnd:
        num.append(float(f_fm[j]))
    pp_f.append(max(num)-min(num))
    num=[]
    rnd=np.random.randint(0,len(m_HJD),int(0.1*len(m_HJD)))
    for j in rnd:
        num.append(float(m_fm[j]))
    pp_m.append(max(num)-min(num)) 

# In[summary]
plt.close('all')
try: plt.close(sp)
except: pass
sp,ax1=plt.subplots(2,2)
ax1[0,0].plot(np.vstack(s_HJD),np.vstack(s_fm)-np.mean(s_fm),'.',label='Rescaled SONG')
ax1[0,0].plot(np.vstack(f_HJD),np.vstack(f_fm)-np.mean(f_fm),'.',label='Rescaled FIES')            
ax1[0,0].plot(np.vstack(m_HJD),np.vstack(m_fm)-np.mean(m_fm),'.',label='Rescaled HERMES')        
ax1[0,0].grid(),ax1[0,0].legend()
ax1[0,0].set_xlabel('HJD [days]-2450000'),ax1[0,0].set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$'),

h_s=ax1[1,0].hist(pp_s,bins=30,density=True,alpha=0.3,label='SONG')
h_f=ax1[1,0].hist(pp_f,bins=30,density=True,alpha=0.3,label='FIES')
h_m=ax1[1,0].hist(pp_m,bins=30,density=True,alpha=0.3,label='HERMES')
ax1[1,0].legend()
ax1[1,0].set_xlabel('Peak-to-peak [km/s]')
ax1[1,0].plot([h_s[1][h_s[0].argmax()]-2*np.mean(error_s),h_s[1][h_s[0].argmax()]+2*np.mean(error_s)],[max(h_s[0]),max(h_s[0])],color='#1f77b4')
ax1[1,0].plot([h_s[1][h_s[0].argmax()]-2*np.mean(error_s),h_s[1][h_s[0].argmax()]+2*np.mean(error_s)],[max(h_s[0]),max(h_s[0])],color='#1f77b4')

ax1[1,0].plot([h_f[1][h_f[0].argmax()]-2*np.mean(error_f),h_f[1][h_f[0].argmax()]+2*np.mean(error_f)],[max(h_f[0]),max(h_f[0])],color='#1f7f0e')
ax1[1,0].plot([h_f[1][h_f[0].argmax()]-2*np.mean(error_f),h_f[1][h_f[0].argmax()]+2*np.mean(error_f)],[max(h_f[0]),max(h_f[0])],color='#1f7f0e')

ax1[1,0].plot([h_m[1][h_m[0].argmax()]-2*np.mean(error_m),h_m[1][h_m[0].argmax()]+2*np.mean(error_m)],[max(h_m[0]),max(h_m[0])],color='#2ca02c')
ax1[1,0].plot([h_m[1][h_m[0].argmax()]-2*np.mean(error_m),h_m[1][h_m[0].argmax()]+2*np.mean(error_m)],[max(h_m[0]),max(h_m[0])],color='#2ca02c')

ax1[0,1].set_title('Random popints')
rn_s=ax1[0,1].hist(STD_s_rn,bins=50,density=True,alpha=0.3,label='SONG')
rn_f=ax1[0,1].hist(STD_f_rn,bins=50,density=True,alpha=0.3,label='FIES')
rn_m=ax1[0,1].hist(STD_m_rn,bins=50,density=True,alpha=0.3,label='HERMES')
ax1[0,1].legend()
from scipy import stats




ax1[1,1].hist(STD_s_cv,bins=20,density=True,alpha=0.3,label='SONG')
ax1[1,1].hist(STD_f_cv,bins=20,density=True,alpha=0.3,label='FIES') 
ax1[1,1].hist(STD_m_cv,bins=20,density=True,alpha=0.3,label='HERMES') 
ax1[1,1].legend()
ax1[1,1].set_xlabel(r'$\sigma_v$ [km/s]')
ax1[1,1].set_title('Convolution')
ax1[0,1].get_shared_x_axes().join(ax1[0,1], ax1[1,1])
ax1[0,1].set_xticklabels([])


# Hay que buscar la moda (máximo del bineado... )

from matplotlib.patches import Rectangle
ax1[0,0].add_patch(Rectangle((data_0, rn_s[1][rn_s[0].argmax()]-np.std(rn_s[1])),data_1-data_0, 2*np.std(rn_s[1]),alpha=0.3,color='#1f77b4'))
ax1[0,0].add_patch(Rectangle((data_0, -rn_s[1][rn_s[0].argmax()]-np.std(rn_s[1])),data_1-data_0, 2*np.std(rn_s[1]),alpha=0.3,color='#1f77b4'))

ax1[0,0].add_patch(Rectangle((data_0, rn_f[1][rn_f[0].argmax()]-np.std(rn_f[1])),data_1-data_0, 2*np.std(rn_f[1]),alpha=0.3,color='#1f7f0e'))
ax1[0,0].add_patch(Rectangle((data_0, -rn_f[1][rn_f[0].argmax()]-np.std(rn_f[1])),data_1-data_0, 2*np.std(rn_f[1]),alpha=0.3,color='#1f7f0e'))

ax1[0,0].add_patch(Rectangle((data_0, rn_m[1][rn_m[0].argmax()]-np.std(rn_m[1])),data_1-data_0, 2*np.std(rn_m[1]),alpha=0.3,color='#2ca02c'))
ax1[0,0].add_patch(Rectangle((data_0, -rn_m[1][rn_m[0].argmax()]-np.std(rn_m[1])),data_1-data_0, 2*np.std(rn_m[1]),alpha=0.3,color='#2ca02c'))


plt.show()




# probar conv para 1 y 10 días, de modo que saquemos cuál es periodo de variabilidad característico 