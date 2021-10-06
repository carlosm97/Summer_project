#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 12:33:18 2021

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

# In[Reading_TESS]
path='../updated_northern_lc/new_lc_ascii'

lc=os.listdir(path)

# HD2905
for j in lc:
    if '419531713' in j:
        print('HD2905 is',j)
        HD2905=j
hdu = open(path+'/'+HD2905)
HJD,mag=[],[]
for j in hdu:
    dat=j.split()
    HJD.append(float(dat[0])+7000.0),mag.append(float(dat[1]))
hdu.close()




SM=reading_summary('../stars/HD2905/4574/00_HD2905_4574.txt')
m_4574_fm,f_4574_fm,s_4574_fm=[],[],[]
m_HJD,f_HJD,s_HJD=[],[],[]
for i in range(len(SM.first_moment)):
    if SM.Instrument[i]=="'MERCATOR'":
        m_4574_fm.append(SM.first_moment[i]),m_HJD.append(SM.HJD[i])
    elif SM.Instrument[i]=="'FIES'":
        f_4574_fm.append(SM.first_moment[i]),f_HJD.append(SM.HJD[i])
    elif SM.Instrument[i]=="'SONG'":
        s_4574_fm.append(SM.first_moment[i]),s_HJD.append(SM.HJD[i])
    else:
        print('Unknown instrument')

Z_m=sorted(zip(np.vstack(m_HJD)-2.45e6,np.vstack(m_4574_fm)-np.mean(m_4574_fm)))
Z_f=sorted(zip(np.vstack(f_HJD)-2.45e6,np.vstack(f_4574_fm)-np.mean(f_4574_fm)))
Z_s=sorted(zip(np.vstack(s_HJD)-2.45e6,np.vstack(s_4574_fm)-np.mean(s_4574_fm)))

m_HJD_4574,m_fm_4574=[Z_m[i][0] for i in range(len(Z_m))],[Z_m[i][1] for i in range(len(Z_m))]
f_HJD_4574,f_fm_4574=[Z_f[i][0] for i in range(len(Z_f))],[Z_f[i][1] for i in range(len(Z_f))]
s_HJD_4574,s_fm_4574=[Z_s[i][0] for i in range(len(Z_s))],[Z_s[i][1] for i in range(len(Z_s))]

plt.close(fig)
fig,ax=plt.subplots(3,sharex=True)
ax[0].plot(HJD,mag,'.',markersize=0.5)
ax[0].grid(),ax[0].set_ylabel('Magnitude')

orf,lf,jd,zm=zero_reading('HD2905','H_alpha_6562')
ax[1].plot(np.vstack(jd)-2.45e6,zm,'.',label=r'H$\alpha$')
ax[1].grid(),ax[1].legend()
ax[1].set_ylabel(r'H$\alpha$')

ax[2].plot(np.vstack(m_HJD_4574),np.vstack(m_fm_4574)-np.mean(m_fm_4574),'.r')#,label='SiIII HERMES')        
ax[2].plot(np.vstack(f_HJD_4574),np.vstack(f_fm_4574)-np.mean(f_fm_4574),'.r')#,label='SiIII FIES')            
ax[2].plot(np.vstack(s_HJD_4574),np.vstack(s_fm_4574)-np.mean(s_fm_4574),'.r',label='SiIII')
ax[2].grid(),ax[2].legend()
ax[2].set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$'),ax[2].set_xlabel('HJD [days]-2450000')

ax[1].set_xlim(8750,9000)

plt.show()

# In[Time_variability]

pp_LC,std_LC=[],[]
for i in range(1000):
    # Selection of the data to be used:
    num=[]
    rnd=np.random.randint(0,len(mag),int(0.1*len(mag)))
    for j in rnd:
        num.append(float(mag[j]))
    pp_LC.append(max(num)-min(num))
    std_LC.append(np.std(num))

fpp, app= plt.subplots(2)
app[0].hist(pp_LC,bins=50)
app[1].hist(std_LC,bins=50)

# In[ODS]

from pandas_ods_reader import read_ods

base_path = "../updated_northern_lc/Updated_notes.ods"
sheet_index = 1
df = read_ods(base_path , sheet_index)
TIC=np.array(df.TIC)
HD=np.array(df.HD)










