#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:08:39 2021

@author: charlie
"""
import numpy as np 
import matplotlib.pyplot as plt
import os
os.chdir('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/pyIACOB')  
from scipy.optimize import curve_fit
from scipy import integrate
import formulas as fmr
from scipy.interpolate import interp1d
from db import *
import spec as sp
from RV import *
from line_study import *
import progressbar

star='HD91316' # B1Iab --> NO HeII
files=findstar(star) # Solo hay espectros FH. 
FL=[]
for f in files: 
    if 'ST' not in f and 'AST' not in f :
        FL.append(f)
print(len(FL))
#%%
# Comrpobamos que la corrección heliocéntrica es correcta. 
fig,axs=plt.subplots()
I=0
bar = progressbar.ProgressBar(max_value=len(FL))

for f in FL: 
    try: 
        SM=sp.spec(f,line=5892.5)
        SM.resamp(0.015,lwl=5887.5,rwl=5897.5)
        axs.plot(SM.wave,SM.flux)
    except:
        print(f)
    I+=1
    bar.update(I)
#5889.950 5895.924 
#axs.set_xlim(4000,7000)
axs.set_ylim(0.,1.2)
axs.grid()
plt.show()
# In[Metals]
plt.close('all')
prueba_linea(star,'FIES',4552.622,4545,4560)

prueba_linea(star,'MERCATOR',4552.622,4545,4560)# telluric=True

prueba_linea(star,'SONG',4552.622,4549,4557)

SiIII_1=4552.622
SiIII_2=4567.840
SiIV=4116.103

HD91316=line_study(star,'FIES',SiIII_1,4549,4557,0.015,plt_mean=True,lim1=0.95,lim2=1.0,telluric=True)
HD91316=line_study(star,'MERCATOR',SiIII_1,4549,4557,0.015,plt_mean=True,lim1=1.05,lim2=1.1,telluric=True)
#HD91316=line_study(star,'SONG',4552.622,4549,4557,0.015,plt_mean=True,lim1=1.05,lim2=1.1,telluric=True,save=False,snr_min=15)

HD91316=line_study(star,'FIES',SiIII_2,4565,4572,0.015,lim2=1.1,plt_mean=True)
HD91316=line_study(star,'MERCATOR',SiIII_2,4565,4573,0.015,lim1=1.1,plt_mean=True,telluric=True)
HD91316=line_study(star,'SONG',SiIII_2,4565,4573,0.015,lim1=1.2,lim2=1.25,plt_mean=True,snr_min=15,telluric=True)

HD91316=line_study(star,'FIES',SiIV,4115,4118.5,0.015,plt_mean=True,telluric=True,automatic_selection=False,ll1=-50,ll2=150,sel='vel')
HD91316=line_study(star,'MERCATOR',SiIV,4115,4118.5,0.015,plt_mean=True,telluric=True,automatic_selection=False,ll1=-50,ll2=150,sel='vel')

OII=4661.6341
HD91316=line_study(star,'FIES',OII,4659,4665,0.015,lim2=1.05,plt_mean=True)
HD91316=line_study(star,'MERCATOR',OII,4659,4665,0.015,lim2=1.05,telluric=True,plt_mean=True)
HD91316=line_study(star,'SONG',OII,4660,4665,0.015,telluric=True,plt_mean=True,snr_min=15,automatic_selection=False,sel='vel',ll1=-50,ll2=150)

OIII=5592.252 # Not
SiII=6347.11 # Not
SiIII_3=5739.734


HD91316=line_study(star,'FIES',SiIII_3,5736,5745,0.015,plt_mean=True,lim1=0.95,lim2=0.96,sig=3,sel='vel')
HD91316=line_study(star,'MERCATOR',SiIII_3,5736,5745,0.015,plt_mean=True,lim1=1,lim2=1.05,telluric=True)
HD91316=line_study(star,'SONG',SiIII_3,5737,5743,0.015,lim1=0.85,lim2=0.9,plt_mean=True,telluric=True,snr_min=15,sig=3,sel='vel')


# In[hydrogen]

H_alpha, H_beta, HeI=6562.80,4861.325,5875.62

prueba_linea(star,'FIES',H_alpha,6540,6585)

HD91316=line_study(star,'FIES',H_beta,4845,4880,0.015,plt_mean=True,target='H_beta',automatic_selection=False,ll1=-500,ll2=650,sel='vel')
HD91316=line_study(star,'MERCATOR',H_beta,4845,4880,0.015,plt_mean=True,target='H_beta',automatic_selection=False,telluric=True,ll1=-500,ll2=500,sel='vel')
HD91316=line_study(star,'SONG',H_beta,4850,4873,0.015,lim1=1.6,lim2=1.8,plt_mean=True,target='H_beta',snr_min=15,telluric=True)

HD91316=line_study(star,'FIES',H_alpha,6550,6575,0.015,plt_mean=True,telluric=True,target='H_alpha')
HD91316=line_study(star,'MERCATOR',H_alpha,6550,6575,0.015,lim1=1.4,lim2=1.6,plt_mean=True,telluric=True,target='H_alpha')
HD91316=line_study(star,'SONG',H_alpha,6556,6570,0.015,plt_mean=True,telluric=True,target='H_alpha',lim1=1.6,lim2=1.7,snr_min=15)


# In[He I]
HeI_1=4387.929
HeI_2=4713.145
HeI_3=5015.678
HeI_4=5875.62
prueba_linea(star,'FIES',HeI_1,4380,4395)
HD91316=line_study(star,'FIES',HeI_1,4385,4392,0.015,plt_mean=True,target='HeI_4387',automatic_selection=False,sel='vel',ll1=-100,ll2=180)
HD91316=line_study(star,'MERCATOR',HeI_1,4385,4392,0.015,plt_mean=True,target='HeI_4387',lim2=1.1)


HD91316=line_study(star,'FIES',HeI_2,4711,4717,0.015,plt_mean=True,target='HeI_4713',automatic_selection=False,ll1=4712,ll2=4716)
HD91316=line_study(star,'MERCATOR',HeI_2,4711,4717,0.015,plt_mean=True,target='HeI_4713',automatic_selection=False,ll1=4712,ll2=4716)
HD91316=line_study(star,'SONG',HeI_2,4711,4717,0.015,plt_mean=True,target='HeI_4713',telluric=True,automatic_selection=False,ll1=4712,ll2=4716,snr_min=15)


HD91316=line_study(star,'FIES',HeI_3,5013.3,5019.5,0.015,lim1=0.9,lim2=0.93,plt_mean=True,target='HeI_5015',telluric=True)
HD91316=line_study(star,'MERCATOR',HeI_3,5013.3,5019.5,0.015,lim2=1.1,plt_mean=True,target='HeI_5015')
HD91316=line_study(star,'SONG',HeI_3,5013.3,5019.5,0.015,lim1=1.3,lim2=1.4,plt_mean=True,target='HeI_5015',snr_min=15,telluric=True)


#HD91316=line_study(star,'FIES',HeI_4,5872,5881,0.015,lim1=0.9,lim2=0.93,plt_mean=True,target='HeI_5015',telluric=True) # A la espera de hablar con Sergio
HD91316=line_study(star,'MERCATOR',HeI_4,5872,5881,0.015,plt_mean=True,target='HeI_5875',telluric=True)
HD91316=line_study(star,'SONG',HeI_4,5872,5881,0.015,lim1=1.2,lim2=1.25,plt_mean=True,target='HeI_5875',telluric=True,snr_min=15)



# In[Strange]

condition='_N_'
ext='_FH'
spectrum=findstar(star)

files=[ins for ins in spectrum if condition in ins]
strange=[]
for s in files:
    SM=sp.spec(s)
    wave,flux=SM.wave,SM.flux
    SM.resamp(0.015,lwl=5874.2,rwl=5874.6)
    if max(SM.flux)>0.98:
        strange.append(s)
        plt.plot(wave,flux)
plt.xlim(5872,5881),plt.ylim(0.5,1.1)
plt.grid()
