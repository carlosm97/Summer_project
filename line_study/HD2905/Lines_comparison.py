#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 20:51:42 2021

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
# In[mean_computation]
runcell(0, '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/HD2905/Lines_comparison.py')
flux_M,date_M,wave_M,v_M=[],[],[],[]
flux_F,date_F,wave_F,v_F=[],[],[],[]
flux_S,date_S,wave_S,v_S=[],[],[],[]
H_alpha=6562.80
WV_M,V_M,F_M,HJD_M,MERCATOR_mean=mean_line('HD2905','H_alpha_6562',H_alpha,'MERCATOR')
WV_F,V_F,F_F,HJD_F,FIES_mean=mean_line('HD2905','H_alpha_6562',H_alpha,'FIES')
WV_S,V_S,F_S,HJD_S,SONG_mean=mean_line('HD2905','H_alpha_6562',H_alpha,'SONG')
fig1,ax1=plt.subplots()
pt=ax1.plot(HJD_M,F_M,'r.')
pt=ax1.plot(HJD_F,F_F,'b.')
pt=ax1.plot(HJD_S,F_S,'g.')

wave,fecha=[],[]
jd1,jd2=2458765,2458795
for i in range(len(HJD_M)):
    if HJD_M[i]>jd1 and HJD_M[i]<jd2:
        flux_M.append(F_M[i]/MERCATOR_mean),date_M.append(HJD_M[i])        
        wave_M.append(WV_M[i]),v_M.append(V_M[i])
        wave.append(WV_M[i]),fecha.append(HJD_M[i])
for i in range(len(HJD_F)):
    if HJD_F[i]>jd1 and HJD_F[i]<jd2:
        flux_F.append(F_F[i]/FIES_mean),date_F.append(HJD_F[i])        
        wave_F.append(WV_F[i]),v_F.append(V_F[i])
        wave.append(WV_F[i]),fecha.append(HJD_F[i])
        
for i in range(len(HJD_S)):
   if HJD_S[i]>jd1 and HJD_S[i]<jd2:
       flux_S.append(F_S[i]/SONG_mean),date_S.append(HJD_S[i])        
       wave_S.append(WV_S[i]),v_S.append(V_S[i])     
       wave.append(WV_S[i]),fecha.append(HJD_S[i]) 
#pt=[plt.plot(WV_M[0],flux[i]) for i in range(len(flux))]
if np.shape(flux_M)!=(0,):
    fecha_M=np.zeros(np.shape(flux_M))
    for i in range(len(fecha_M[0,:])):
        for j in range(len(fecha_M[:,0])): 
            fecha_M[j,i]=date_M[j]
if np.shape(flux_S)!=(0,):
    fecha_S=np.zeros(np.shape(flux_S))
    for i in range(len(fecha_S[0,:])):
        for j in range(len(fecha_S[:,0])): 
            fecha_S[j,i]=date_S[j]
if np.shape(flux_F)!=(0,):
    fecha_F=np.zeros(np.shape(flux_F))
    for i in range(len(fecha_F[0,:])):
        for j in range(len(fecha_F[:,0])): 
            fecha_F[j,i]=date_F[j]


plt.close('all')
cm = plt.cm.get_cmap('RdYlBu')
try: pt = plt.scatter(wave_S,fecha_S,c=flux_S,s=0.5,cmap=cm,vmax=1.05)
except: pass
try: pt = plt.scatter(wave_M,fecha_M,c=flux_M,s=0.5,cmap=cm,vmax=1.05)
except: pass
try: pt = plt.scatter(wave_F,fecha_F,c=flux_F,s=0.5,cmap=cm,vmax=1.05)
except: pass
plt.show()
plt.colorbar()
plt.title(r'Residues, $\lambda_0=$'+str(6562.80))
plt.ylabel('HJD [days]')
plt.xlabel(r'$\lambda(\AA)$')
plt.xlim(wave_S[0][0],wave_S[0][-1])


# In[GIF_H_alpha]
'''
# Long panel
plt.subplot(grid[0, 0])
plt.subplot(grid[0, 1])
plt.subplot(grid[1, 0])
plt.subplot(grid[1, 1])
plt.subplot(grid[2, :])
'''
jd1,jd2=2458765,2458795

plt.close('all')
runcell(0, '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/HD2905/Lines_comparison.py')
# Añadir los datos del resto de espectrógrafos+trabajar con el intervalo temporal que me interesa

gif('HD2905', 'H_alpha_6562', date_S, F_S, v_S, SONG_mean, jd1, jd2)


# In[Resume]
try: plt.close(fig)
except: pass
fig,(ax1,ax2,ax3)=plt.subplots(3,sharex=True) # Share x axis
# Reference lines:
ax1.plot(np.linspace(-15,15,100),np.linspace(-15,15,100),'--',c='grey')
ax2.plot(np.linspace(-15,15,100),np.linspace(-15,15,100),'--',c='grey')
ax3.plot(np.linspace(-15,15,100),np.linspace(-15,15,100),'--',c='grey')

# First subplot
m_SiIII_OII1,f_SiIII_OII1,s_SiIII_OII1=fm_fm('4590')
m_SiIII_OII2,f_SiIII_OII2,s_SiIII_OII2=fm_fm('4661')
m_SiIII_SiIV,f_SiIII_SiIV,s_SiIII_SiIV=fm_fm('4116')
m_SiIII_SiIII2,f_SiIII_SiIII2,s_SiIII_SiIII2=fm_fm('4567')
m_SiIII_SiIII3,f_SiIII_SiIII3,s_SiIII_SiIII3=fm_fm('4574')


ax1.plot(m_SiIII_SiIII2[:,0]-np.mean(m_SiIII_SiIII2[:,0]),m_SiIII_SiIII2[:,1]-np.mean(m_SiIII_SiIII2[:,1]),'.k')
ax1.plot(f_SiIII_SiIII2[:,0]-np.mean(f_SiIII_SiIII2[:,0]),f_SiIII_SiIII2[:,1]-np.mean(f_SiIII_SiIII2[:,1]),'xk',label=r'Si III (4567.840 $\AA$)')
ax1.plot(m_SiIII_SiIII3[:,0]-np.mean(m_SiIII_SiIII3[:,0]),m_SiIII_SiIII3[:,1]-np.mean(m_SiIII_SiIII3[:,1]),'.r')
ax1.plot(f_SiIII_SiIII3[:,0]-np.mean(f_SiIII_SiIII3[:,0]),f_SiIII_SiIII3[:,1]-np.mean(f_SiIII_SiIII3[:,1]),'xr',label=r'Si III (4574.757 $\AA$)')
ax1.plot(m_SiIII_SiIV[:,0]-np.mean(m_SiIII_SiIV[:,0]),m_SiIII_SiIV[:,1]-np.mean(m_SiIII_SiIV[:,1]),'.b')
ax1.plot(f_SiIII_SiIV[:,0]-np.mean(f_SiIII_SiIV[:,0]),f_SiIII_SiIV[:,1]-np.mean(f_SiIII_SiIV[:,1]),'xb',label=r'Si IV (4116.103 $\AA$)')
ax1.plot(m_SiIII_OII2[:,0]-np.mean(m_SiIII_OII2[:,0]),m_SiIII_OII2[:,1]-np.mean(m_SiIII_OII2[:,1]),'.g')
ax1.plot(f_SiIII_OII2[:,0]-np.mean(f_SiIII_OII2[:,0]),f_SiIII_OII2[:,1]-np.mean(f_SiIII_OII2[:,1]),'xg',label=r'O II (4661.6324 $\AA$)')
ax1.plot(m_SiIII_OII1[:,0]-np.mean(m_SiIII_OII1[:,0]),m_SiIII_OII1[:,1]-np.mean(m_SiIII_OII1[:,1]),'.',c='orange')
ax1.plot(f_SiIII_OII1[:,0]-np.mean(f_SiIII_OII1[:,0]),f_SiIII_OII1[:,1]-np.mean(f_SiIII_OII1[:,1]),'x',c='orange',label=r'O II (4590.974 $\AA$)')


# Second subplot
m_SiIII_HeI_5875,f_SiIII_HeI_5875,s_SiIII_HeI_5875=fm_fm('HeI_5875')
m_SiIII_HeI_4387,f_SiIII_HeI_4387,s_SiIII_HeI_4387=fm_fm('HeI_4387')
m_SiIII_HeI_5015,f_SiIII_HeI_5015,s_SiIII_HeI_5015=fm_fm('HeI_5015')

ax2.plot(m_SiIII_HeI_4387[:,0]-np.mean(m_SiIII_HeI_4387[:,0]),m_SiIII_HeI_4387[:,1]-np.mean(m_SiIII_HeI_4387[:,1]),'.k')
ax2.plot(f_SiIII_HeI_4387[:,0]-np.mean(f_SiIII_HeI_4387[:,0]),f_SiIII_HeI_4387[:,1]-np.mean(f_SiIII_HeI_4387[:,1]),'xk',label=r'He I (4387.929 $\AA$)')
ax2.plot(m_SiIII_HeI_5015[:,0]-np.mean(m_SiIII_HeI_5015[:,0]),m_SiIII_HeI_5015[:,1]-np.mean(m_SiIII_HeI_5015[:,1]),'.r')
ax2.plot(f_SiIII_HeI_5015[:,0]-np.mean(f_SiIII_HeI_5015[:,0]),f_SiIII_HeI_5015[:,1]-np.mean(f_SiIII_HeI_5015[:,1]),'xr',label=r'He I (5015.678$\AA$)')
ax2.plot(m_SiIII_HeI_5875[:,0]-np.mean(m_SiIII_HeI_5875[:,0]),m_SiIII_HeI_5875[:,1]-np.mean(m_SiIII_HeI_5875[:,1]),'.b')
ax2.plot(f_SiIII_HeI_5875[:,0]-np.mean(f_SiIII_HeI_5875[:,0]),f_SiIII_HeI_5875[:,1]-np.mean(f_SiIII_HeI_5875[:,1]),'xb',label=r'He I (5875.62 $\AA$)')

# Third subplot
m_SiIII_Halpha,f_SiIII_Halpha,s_SiIII_Halpha=fm_fm('H_alpha_6562')
m_SiIII_Hbeta,f_SiIII_Hbeta,s_SiIII_Hbeta=fm_fm('H_beta_4861')

ax3.plot(m_SiIII_Halpha[:,0]-np.mean(m_SiIII_Halpha[:,0]),-0.5*(m_SiIII_Halpha[:,1]-np.mean(m_SiIII_Halpha[:,1])),'.b')
ax3.plot(f_SiIII_Halpha[:,0]-np.mean(f_SiIII_Halpha[:,0]),-0.5*(f_SiIII_Halpha[:,1]-np.mean(f_SiIII_Halpha[:,1])),'xb',label=r'-0.5$\cdot$H $\alpha$ (6562.8 $\AA$)')
ax3.plot(m_SiIII_Hbeta[:,0]-np.mean(m_SiIII_Hbeta[:,0]),m_SiIII_Hbeta[:,1]-np.mean(m_SiIII_Hbeta[:,1]),'.r')
ax3.plot(f_SiIII_Hbeta[:,0]-np.mean(f_SiIII_Hbeta[:,0]),f_SiIII_Hbeta[:,1]-np.mean(f_SiIII_Hbeta[:,1]),'xr',label=r'H $\beta$ (4861.32 $\AA$)')

# General and format things 
ax1.grid(),ax2.grid(),ax3.grid()
ax1.legend(),ax2.legend(),ax3.legend()
ax3.set_xlabel(r'$\langle v\rangle-\langle v_0\rangle$ Si III (4552.622) [km/s]')
ax1.set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ [km/s]'),ax2.set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ [km/s]'),ax3.set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ [km/s]')

try: plt.close(fig2)
except: pass

# ZONES TO ZOOM IN 
l1,l2,l3,l4=6640,6650,7575,7620
fig2,axs=plt.subplots(2,3,sharey='row', sharex='col')
axs[1,1].set_xlabel('HJD [days]-2450000')
folders='4552','H_alpha_6562','H_beta_4861','HeI_5875','4661'
star='HD2905'

OF_SiIII,LF_SiIII,HJD_SiIII,ZM_SiIII=zero_reading(star,folders[0])
OF_Halpha,LF_Halpha,HJD_Halpha,ZM_Halpha=zero_reading(star,folders[1])
OF_Hbeta,LF_Hbeta,HJD_Hbeta,ZM_Hbeta=zero_reading(star,folders[2])
OF_HeI,LF_HeI,HJD_HeI,ZM_HeI=zero_reading(star,folders[3])
OF_OII,LF_OII,HJD_OII,ZM_OII=zero_reading(star,folders[4])

axs[0,0].plot(HJD_SiIII-2.45e6,ZM_SiIII-np.mean(ZM_SiIII),'s',label='Si III')
axs[0,0].plot(HJD_OII-2.45e6,ZM_OII-np.mean(ZM_OII),'^',label='O II')
axs[0,0].plot(HJD_Halpha-2.45e6,ZM_Halpha-np.mean(ZM_Halpha),'.',label=r'H$\alpha$')
axs[0,0].plot(HJD_Hbeta-2.45e6,ZM_Hbeta-np.mean(ZM_Hbeta),'.',label=r'H$\beta$')
axs[0,0].plot(HJD_HeI-2.45e6,ZM_HeI-np.mean(ZM_HeI),'.',label='He I')
axs[0,0].grid(),axs[0,0].legend(),axs[0,0].set_ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
axs[0,0].axvline(x=l1,c='k'),axs[0,0].axvline(x=l2,c='k'),axs[0,0].axvline(x=l3,c='k'),axs[0,0].axvline(x=l4,c='k')

axs[0,1].plot(HJD_SiIII-2.45e6,ZM_SiIII-np.mean(ZM_SiIII),'s',label='Si III')
axs[0,1].plot(HJD_OII-2.45e6,ZM_OII-np.mean(ZM_OII),'^',label='O II')
axs[0,1].plot(HJD_Halpha-2.45e6,ZM_Halpha-np.mean(ZM_Halpha),'.',label=r'H$\alpha$')
axs[0,1].plot(HJD_Hbeta-2.45e6,ZM_Hbeta-np.mean(ZM_Hbeta),'.',label=r'H$\beta$')
axs[0,1].plot(HJD_HeI-2.45e6,ZM_HeI-np.mean(ZM_HeI),'.',label='He I')
axs[0,1].grid(),axs[0,1].legend()#,axs[0,1].set_ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
axs[0,1].set_xlim(l1,l2)#,axs[0,1].set_ylim(-75,25)


axs[0,2].plot(HJD_SiIII-2.45e6,ZM_SiIII-np.mean(ZM_SiIII),'s',label='Si III')
axs[0,2].plot(HJD_OII-2.45e6,ZM_OII-np.mean(ZM_OII),'^',label='O II')
axs[0,2].plot(HJD_Halpha-2.45e6,ZM_Halpha-np.mean(ZM_Halpha),'.',label=r'H$\alpha$')
axs[0,2].plot(HJD_Hbeta-2.45e6,ZM_Hbeta-np.mean(ZM_Hbeta),'.',label=r'H$\beta$')
axs[0,2].plot(HJD_HeI-2.45e6,ZM_HeI-np.mean(ZM_HeI),'.',label='He I')
axs[0,2].grid(),axs[0,2].legend()#,axs[0,2].set_xlabel('HJD [days]-2450000')#,axs[0,2].set_ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
axs[0,2].set_xlim(l3,l4)#,axs[0,2].set_ylim(-50,50)


SiIII=reading_summary('../stars/HD2905/4552/00_HD2905_4552.txt')
SiIV=reading_summary('../stars/HD2905/4116/00_HD2905_4116.txt')
OII=reading_summary('../stars/HD2905/4661/00_HD2905_4661.txt')

axs[1,0].plot(np.vstack(SiIII.HJD)-2.45e6,SiIII.first_moment,'s',label='Si III')
axs[1,0].plot(np.vstack(SiIV.HJD)-2.45e6,SiIV.first_moment,'.',label='Si IV')
axs[1,0].plot(np.vstack(OII.HJD)-2.45e6,OII.first_moment,'.',label='O II')
axs[1,0].legend(),axs[1,0].grid(),axs[1,0].set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ [km/s]')
axs[1,0].axvline(x=l1,c='k'),axs[1,0].axvline(x=l2,c='k'),axs[1,0].axvline(x=l3,c='k'),axs[1,0].axvline(x=l4,c='k')

axs[1,1].plot(np.vstack(SiIII.HJD)-2.45e6,SiIII.first_moment,'s',label='Si III')
axs[1,1].plot(np.vstack(SiIV.HJD)-2.45e6,SiIV.first_moment,'.',label='Si IV')
axs[1,1].plot(np.vstack(OII.HJD)-2.45e6,OII.first_moment,'.',label='O II')
axs[1,1].legend(),axs[1,1].grid()

axs[1,2].plot(np.vstack(SiIII.HJD)-2.45e6,SiIII.first_moment,'s',label='Si III')
axs[1,2].plot(np.vstack(SiIV.HJD)-2.45e6,SiIV.first_moment,'.',label='Si IV')
axs[1,2].plot(np.vstack(OII.HJD)-2.45e6,OII.first_moment,'.',label='O II')
axs[1,2].legend(),axs[1,2].grid()

# In[Selection_FIES]

# Buscamos qué sucede con Si IV en FIES. 

path='../stars/HD2905/4116/'

SM=reading_summary(path+'00_HD2905_4116.txt')
minor,maior=[],[]
file_minor,file_maior=[],[]
for i in range(len(SM.Instrument)):
    if SM.Instrument[i]=="'FIES'":
        if SM.first_moment[i]<-5:
            minor.append(SM.first_moment[i]),file_minor.append(SM.line_file[i][1:-1])
        else: 
            maior.append(SM.first_moment[i]),file_maior.append(SM.line_file[i][1:-1])
fig,axis=plt.subplots()
plt.plot()

axis.plot(minor,'.'),axis.plot(maior,'.')
axis.grid()









