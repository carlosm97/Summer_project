#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 15:55:48 2021

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

star='HD34078' # O9.5V
files=findstar(star) # Solo hay espectros FH. 
#%%
# Comrpobamos que la corrección heliocéntrica es correcta. 
fig,axs=plt.subplots()
for f in files: 
    SM=sp.spec(f)
    axs.plot(SM.wave,SM.flux)
axs.set_xlim(5887.5,5897.5)
#axs.set_xlim(4000,7000)
axs.set_ylim(0.,1.2)
axs.grid()
plt.show()
# In[looking_for_lines]
# Probamos las líneas clave: 
SiIII=4552.622
OIII=5592.252
SiII=6347.11 # No presence
OII_1=4661.6324

SiIV=4116.103
SiIII_2=4567.840
SiIII_3=4574.757

NV=4603.74
prueba_linea(star,'MERCATOR',NV,4600,4610)
prueba_linea(star,'FIES',NV,4600,4610)


prueba_linea(star,'FIES',SiIV,4115,4119) # Hay 2 espectros que se van, pero parecen bien corregidos en el resto... 

prueba_linea(star,'MERCATOR',SiIII,4549,4556)
prueba_linea(star,'FIES',SiIII,4549,4556)

prueba_linea(star,'MERCATOR',SiII,6340,6360)
prueba_linea(star,'FIES',SiII,6340,6360)

prueba_linea(star,'MERCATOR',OIII,5580,5600)
prueba_linea(star,'FIES',OIII,5580,5600)

prueba_linea(star,'MERCATOR',OIII,5580,5600)
prueba_linea(star,'FIES',OIII,5580,5600)

prueba_linea(star,'MERCATOR',OII_1,4550,4670)
prueba_linea(star,'FIES',OII_1,4550,4670)

H_alpha, H_beta, HeI=6562.80,4861.325,5875.62

HeII=5411.52

prueba_linea(star,'MERCATOR',H_alpha,6555,6575)
prueba_linea(star,'FIES',H_alpha,6540,6590)

prueba_linea(star,'MERCATOR',H_beta,4840,4890)
prueba_linea(star,'FIES',H_beta,4840,4890)

HeII_2=4859.319

prueba_linea(star,'MERCATOR',HeII_2,4830,4900)
prueba_linea(star,'FIES',HeII_2,4830,4900)

SiIII_2=4567.840
prueba_linea(star,'FIES',SiIII_2,4560,4575) # Preguntar a Sergio por línea verde 
# In[Metals]
OIII=5592.252
HD34078=line_study(star,'FIES',OIII,5590,5596,0.015,plt_mean=True,automatic_selection=False,ll1=-30,ll2=120,sel='vel')
HD34078=line_study(star,'MERCATOR',OIII,5590,5596,0.015,plt_mean=True,automatic_selection=False,ll1=-30,ll2=120,sel='vel')

OII_1=4661.6324
HD34078=line_study(star,'FIES',OII_1,4661.3,4663.5,0.015,plt_mean=True,automatic_selection=False,ll1=10,ll2=100,sel='vel')
HD34078=line_study(star,'MERCATOR',OII_1,4661.3,4663.5,0.015,plt_mean=True,automatic_selection=False,ll1=20,ll2=90,sel='vel')

SiIII=4552.622
HD34078=line_study(star,'FIES',SiIII,4552,4555,0.015,plt_mean=True,automatic_selection=False,ll1=5,ll2=120,sel='vel')
HD34078=line_study(star,'MERCATOR',SiIII,4552,4555,0.015,plt_mean=True,automatic_selection=False,ll1=15,ll2=100,sel='vel')
SiIV=4116.103
HD34078=line_study(star,'FIES',SiIV,4115,4118.5,0.015,plt_mean=True,automatic_selection=False,ll1=4116,ll2=4117.7)
HD34078=line_study(star,'MERCATOR',SiIV,4115,4118.5,0.015,plt_mean=True,automatic_selection=False,ll1=4116,ll2=4117.7)
'''
NV=4603.74
HD34078=line_study(star,'FIES',NV,4600,4606,0.015,plt_mean=True,automatic_selection=False,ll1=4602,ll2=4604,telluric=True)
'''
HD34078=line_study(star,'MERCATOR',SiIII_2,4566,4570.5,0.015,plt_mean=True,automatic_selection=False,ll1=4567.7,ll2=4569.5)
HD34078=line_study(star,'FIES',SiIII_2,4566,4570.5,0.015,plt_mean=True,automatic_selection=False,ll1=4567.7,ll2=4569.5)


# In[Hydrogen]
H_alpha, H_beta=6562.80,4861.325 
HD34078=line_study(star,'FIES',H_beta,4835,4890,0.015,plt_mean=True,automatic_selection=False,ll1=-1000,ll2=1000,sel='vel',target='H_beta_4861')
HD34078=line_study(star,'MERCATOR',H_beta,4835,4890,0.015,plt_mean=True,automatic_selection=False,ll1=-1000,ll2=1000,sel='vel',target='H_beta_4861',telluric=True)

HD34078=line_study(star,'FIES',H_alpha,6530,6600,0.015,plt_mean=True,automatic_selection=False,ll1=-1000,ll2=1000,sel='vel',telluric=True,target='H_alpha')
HD34078=line_study(star,'MERCATOR',H_alpha,6530,6600,0.015,plt_mean=True,automatic_selection=False,ll1=-1000,ll2=1000,sel='vel',telluric=True,target='H_alpha')

# In[HeII]
HeII=5411.52
HD34078=line_study(star,'FIES',HeII,5400,5425,0.015,plt_mean=True,automatic_selection=False,ll1=5407.5,ll2=5417.5,target='HeII')
HD34078=line_study(star,'MERCATOR',HeII,5400,5425,0.015,plt_mean=True,automatic_selection=False,ll1=5407.5,ll2=5417.5,target='HeII',telluric=True)

HeII_2=4859.319

HD34078=line_study(star,'FIES',HeII_2,4830,4890,0.015,plt_mean=True,automatic_selection=False,ll1=4840,ll2=4880,target='HeII_2',telluric=True)
HD34078=line_study(star,'MERCATOR',HeII_2,4830,4890,0.015,plt_mean=True,automatic_selection=False,ll1=4840,ll2=4880,target='HeII_2',telluric=True)

# Por qué tan malo? Quizás por el centro movido?
# Inventamos una línea pseudo-betra medio corregida de RV. 
#H_pseudoBeta=4861.325*(1+60/299792.458)
#HD34078=line_study(star,'FIES',H_pseudoBeta,4835,4890,0.015,plt_mean=True,automatic_selection=False,ll1=-1000,ll2=900,sel='vel',save=False,target='H_pseudo')

#ll1=6550,ll2=6580

# In[HeI]
HeI_1=4387.929
HeI_2=4713.145
HeI_3=5015.678
HeI_4=5875.62
#(star,'FIES',HeI_1,4380,4395)
HD34078=line_study(star,'FIES',HeI_1,4385,4392,0.015,plt_mean=True,target='HeI_4387',automatic_selection=False,sel='vel',ll1=-100,ll2=180)
HD34078=line_study(star,'MERCATOR',HeI_1,4385,4392,0.015,plt_mean=True,target='HeI_4387',automatic_selection=False,sel='vel',ll1=-100,ll2=180)

HD34078=line_study(star,'FIES',HeI_2,4711,4717,0.015,plt_mean=True,target='HeI_4713',automatic_selection=False,sel='vel',ll1=0,ll2=125)
HD34078=line_study(star,'MERCATOR',HeI_2,4711,4717,0.015,plt_mean=True,target='HeI_4713',automatic_selection=False,sel='vel',ll1=0,ll2=125)


HD34078=line_study(star,'FIES',HeI_3,5013.3,5019.5,0.015,plt_mean=True,target='HeI_5015',automatic_selection=False,sel='vel',ll1=-10,ll2=125)
HD34078=line_study(star,'MERCATOR',HeI_3,5013.3,5019.5,0.015,lim2=1.1,plt_mean=True,target='HeI_5015',automatic_selection=False,sel='vel',ll1=-10,ll2=125)


#HD91316=line_study(star,'FIES',HeI_4,5872,5881,0.015,lim1=0.9,lim2=0.93,plt_mean=True,target='HeI_5015',telluric=True) # A la espera de hablar con Sergio. Problema similar a con HD91316
HD34078=line_study(star,'MERCATOR',HeI_4,5872,5881,0.015,plt_mean=True,target='HeI_5875',telluric=True,automatic_selection=False,sel='vel',ll1=-50,ll2=150)
#HD34078=line_study(star,'SONG',HeI_4,5872,5881,0.015,lim1=1.2,lim2=1.25,plt_mean=True,target='HeI_5875',telluric=True,save=False,snr_min=15,automatic_selection=False,sel='vel',ll1=0,ll2=125,save=False)





# In[Zero_moment]
folders='4552','H_alpha','H_beta_4861','HeII','4661','HeII_2'
lines=4552.622,6562.80,4861.325,5411.52,4661.6324,4859.319
bar = progressbar.ProgressBar(max_value=len(lines))
I=0
for i in range(len(lines)):
    try: os.remove('../stars/'+star+'/'+folders[i]+'/'+'00_zero_moment.txt')
    except: pass
    SM=reading_summary('../stars/'+star+'/'+folders[i]+'/00_'+star+'_'+folders[i]+'.txt')
    files=os.listdir('../stars/'+star+'/'+folders[i])
    f=open('../stars/'+star+'/'+folders[i]+'/'+'00_zero_moment.txt', 'a')
    for fil in range(len(SM.line_file)):
        wv,v,fl=reading_line_file('../stars/'+star+'/'+folders[i]+'/'+SM.line_file[fil][1:-1])
        info=SM.Original_file[fil][2:-1],SM.line_file[fil][1:-1],SM.HJD[fil],zero_moment(wv,fl,lines[i])
        f.writelines(str(info))
        f.writelines('\n')
    I+=1
    bar.update(I)
    f.close()
# In[Comparison]
from line_study import *
cM,eM,cF,eF,cS,eS=fm_fm('4661',line1='5592',star='HD34078',error=True)
plt.close('all')
fig,axs=plt.subplots()
axs.errorbar(cF[:,0]-np.mean(cF[:,0]),cF[:,1]-np.mean(cF[:,1]),xerr=eF[:,0],yerr=eF[:,1],fmt='.',alpha=0.3)
axs.grid(),axs.set_xlabel(r'$\langle v\rangle-\langle v_0\rangle$ O III (5592)'),axs.set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ O II (4661)')
axs.plot(cF[:,0]-np.mean(cF[:,0]),cF[:,1]-np.mean(cF[:,1]),'rx')



cM,eM,cF,eF,cS,eS=fm_fm('4116',line1='4552',star='HD34078',error=True)
fig1,axs1=plt.subplots()
axs1.errorbar(cF[:,0]-np.mean(cF[:,0]),cF[:,1]-np.mean(cF[:,1]),xerr=eF[:,0],yerr=eF[:,1],fmt='.',alpha=0.3)
axs1.grid(),axs1.set_xlabel(r'$\langle v\rangle-\langle v_0\rangle$ Si III (4552)'),axs1.set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ Si IV (4116)')
axs1.plot(cF[:,0]-np.mean(cF[:,0]),cF[:,1]-np.mean(cF[:,1]),'rx')
plt.show()

# In[GIF]
fig,axs=plt.subplots()
axs.plot(photometry(star)[0],photometry(star)[1],'.')
axs.grid(),plt.show()

folder='H_alpha'
plt.figure()
plt.plot(zero_reading(star,folder=folder)[2],zero_reading(star,folder=folder)[3],'.')

plt.close('all')


WV_F,V_F,F_F,HJD_F,FIES_mean=mean_line(star,'H_alpha',H_alpha,'FIES')
flux_F,date_F,wave_F,v_F=[],[],[],[]
wave,fecha=[],[]
jd1,jd2=2458835,2458845
for i in range(len(HJD_F)):
    if HJD_F[i]>jd1 and HJD_F[i]<jd2:
        flux_F.append(F_F[i]/FIES_mean),date_F.append(HJD_F[i])        
        wave_F.append(WV_F[i]),v_F.append(V_F[i])
        wave.append(WV_F[i]),fecha.append(HJD_F[i])
        
gif(star,folder,date_F,F_F,v_F,FIES_mean,jd1,jd2)












