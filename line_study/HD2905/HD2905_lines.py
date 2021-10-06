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

star='HD2905' #B1Ia C
# In[SiIII1]
SiIII=4552.622

HD2905=line_study('HD2905','FIES',SiIII,4549,4556,0.015,lim1=1.1,lim2=1.2,plt_mean=True)      
HD2905=line_study('HD2905','MERCATOR',SiIII,4549,4556,0.015,plt_mean=True,telluric=True)    

#prueba_linea('HD2905','SONG',SiIII,4549,4556)
#HD2905=line_study('HD2905','SONG',SiIII,4549,4556,0.015,snr_min=10,MC=50)    
# Está muy al extremo

# In[SiIII2]

SiIII2=5739.734

HD2905=line_study('HD2905','MERCATOR',SiIII2,5732.5,5747.5,0.015,plt_mean=True,telluric=True)#,lim1=0.85,lim2=1,plt_mean=True)    
HD2905=line_study('HD2905','FIES',SiIII2,5732.5,5747.5,0.015,lim1=1.1,lim2=1.2,plt_mean=True,automatic_selection=False,ll1=-200,ll2=200,sel='vel',telluric=True)
#prueba_linea('HD2905','SONG',SiIII2,5732.5,5747) # Al límite

# In[H_alpha]

H_alpha=6562.80

#prueba_linea('HD2905','FIES',H_alpha,6550,6580)
#prueba_linea('HD2905','MERCATOR',H_alpha,6550,6580)
#HD2905=line_study('HD2905','FIES',HeI,5865,5885,0.015,lim1=0.85,lim2=1,plt_mean=True,target='HeI_5875')    
#HD2905=line_study('HD2905','MERCATOR',HeI,5865,5885,0.015,lim1=0.8,lim2=1,plt_mean=True,target='HeI_5875')
HD2905=line_study('HD2905','SONG',H_alpha,6553,6574,0.015,target='H_alpha_6562',plt_mean=True,snr_min=8,telluric=True)
HD2905=line_study('HD2905','FIES',H_alpha,6545,6580,0.015,lim1=2,lim2=2.5,target='H_alpha_6562',plt_mean=True,telluric=True)#,lim1=0.85,lim2=1,plt_mean=True)    
HD2905=line_study('HD2905','MERCATOR',H_alpha,6550,6574,0.015,lim1=1.5,lim2=2.,target='H_alpha_6562',plt_mean=True,telluric=True)#,lim1=0.85,lim2=1,plt_mean=True)

# In[H_beta]

H_beta=4861.325 
#prueba_linea('HD2905','FIES',H_beta,4850,4870)
#prueba_linea('HD2905','MERCATOR',H_beta,4850,4870)
HD2905=line_study('HD2905','FIES',H_beta,4850,4870,0.015,lim1=1.35,lim2=1.45,target='H_beta_4861',plt_mean=True)
HD2905=line_study('HD2905','MERCATOR',H_beta,4850,4870,0.015,lim1=1.2,lim2=1.4,target='H_beta_4861',plt_mean=True)#,lim1=0.85,lim2=1,plt_mean=True)
HD2905=line_study('HD2905','SONG',H_beta,4850,4870,0.015,lim1=1.2,lim2=1.4,target='H_beta_4861',plt_mean=True,snr_min=15)

# In[SiIV]

SiIV=4116.103

#prueba_linea('HD2905','MERCATOR',SiIV,4114,4118)
#prueba_linea('HD2905','FIES',SiIV,4114,4118)
HD2905=line_study('HD2905','MERCATOR',SiIV,4114,4118,0.015,plt_mean=True,lim1=0.65,lim2=0.7,telluric=True)
HD2905=line_study('HD2905','FIES',SiIV,4114,4118,0.015,plt_mean=True,lim1=0.6,lim2=0.85)
# Para SONG está fuera del rango 

# In[SiII_3_2]

SiIII_2=4567.840

#prueba_linea('HD2905','MERCATOR',SiIII_2,4563,4572)
#prueba_linea('HD2905','FIES',SiIII_2,4563,4572)
HD2905=line_study('HD2905','MERCATOR',SiIII_2,4563,4572,0.015,plt_mean=True,lim1=1.1,lim2=1.25,telluric=True)
HD2905=line_study('HD2905','FIES',SiIII_2,4563,4572,0.015,lim1=1.1,lim2=1.25,plt_mean=True)
HD2905=line_study('HD2905','SONG',SiIII_2,4563,4572,0.015,lim1=1.05,lim2=1.15,plt_mean=True,snr_min=15,c1=200,c2=-150,telluric=True)
# In[SiII_3_3]

SiIII_3=4574.757

#prueba_linea('HD2905','MERCATOR',SiIII_3,4571,4578)
#prueba_linea('HD2905','FIES',SiIII_3,4571,4578)
#prueba_linea('HD2905','SONG',SiIII_3,4571,4578)

HD2905=line_study('HD2905','MERCATOR',SiIII_3,4571,4578,0.015,lim1=1.1,lim2=1.2,plt_mean=True,telluric=True)
HD2905=line_study('HD2905','FIES',SiIII_3,4571,4578,0.015,lim1=1.1,lim2=1.15,plt_mean=True)
HD2905=line_study('HD2905','SONG',SiIII_3,4571,4578,0.015,lim1=1.1,lim2=1.15,plt_mean=True,snr_min=15,telluric=True)

# In[OII_1]
OII_1=4661.6324

#prueba_linea('HD2905','MERCATOR',OII_1,4658,4665)
#prueba_linea('HD2905','FIES',OII_1,4658,4665)
HD2905=line_study('HD2905','MERCATOR',OII_1,4658,4665,0.015,lim1=1.05,plt_mean=True,telluric=True)
HD2905=line_study('HD2905','FIES',OII_1,4658,4665,0.015,plt_mean=True,telluric=True)
HD2905=line_study('HD2905','SONG',OII_1,4658,4665,0.015,plt_mean=True,snr_min=15,lim1=1.2,lim2=1.25,telluric=True)
# In[OII_2]

OII_2=4590.974

#prueba_linea('HD2905','FIES',OII_2,4588,4594)
#prueba_linea('HD2905','MERCATOR',OII_2,4588,4594)
#prueba_linea('HD2905','SONG',OII_2,4588,4594)
HD2905=line_study('HD2905','MERCATOR',OII_2,4588,4594,0.015,lim1=1.05,lim2=1.15,plt_mean=True,telluric=True)
HD2905=line_study('HD2905','FIES',OII_2,4588,4594,0.015,lim1=1.05,lim2=1.15,plt_mean=True,telluric=True)
HD2905=line_study('HD2905','SONG',OII_2,4588,4594,0.015,plt_mean=True,snr_min=15,lim1=1.2,lim2=1.25,telluric=True)


# In[HeI]
HeI=5875.62

HD2905=line_study('HD2905','FIES',HeI,5865,5885,0.015,plt_mean=True,target='HeI_5875',telluric=True)    
HD2905=line_study('HD2905','MERCATOR',HeI,5865,5885,0.015,plt_mean=True,target='HeI_5875',telluric=True)
HD2905=line_study('HD2905','SONG',HeI,5865,5885,0.015,lim1=1.3,lim2=1.5,plt_mean=True,target='HeI_5875',snr_min=10,telluric=True)

He_I1,He_I2,He_I3,He_I4=4387.929,5015.678,4471.489,4713.139
#prueba_linea('HD2905','FIES',He_I1,4384,4391.5)
#prueba_linea('HD2905','MERCATOR',He_I1,4384,4391.5)
#prueba_linea('HD2905','SONG',He_I1,4384,4391.5)

HD2905=line_study('HD2905','MERCATOR',He_I1,4384,4391.5,0.015,lim1=1.1,plt_mean=True,target='HeI_4387',telluric=True)
HD2905=line_study('HD2905','FIES',He_I1,4384,4391.5,0.015,lim1=1.25,lim2=1.35,plt_mean=True,target='HeI_4387',telluric=True)

HD2905=line_study('HD2905','FIES',He_I2,5007.5,5022.5,0.015,lim1=1.05,plt_mean=True,target='HeI_5015',automatic_selection=False,ll1=-200,ll2=180,sel='vel',telluric=True)
HD2905=line_study('HD2905','MERCATOR',He_I2,5007.5,5022.5,0.015,lim1=1.05,plt_mean=True,target='HeI_5015',automatic_selection=False,ll1=-200,ll2=180,sel='vel',telluric=True)
HD2905=line_study('HD2905','SONG',He_I2,5007.5,5022.5,0.015,lim1=1.2,lim2=1.4,plt_mean=True,target='HeI_5015',snr_min=15,telluric=True)


HD2905=line_study('HD2905','SONG',He_I3,4468,4475,0.015,lim1=1.2,lim2=1.4,plt_mean=True,target='HeI_4471',snr_min=15,telluric=True)
HD2905=line_study('HD2905','MERCATOR',He_I3,4468,4475,0.015,lim1=1.2,lim2=1.4,plt_mean=True,target='HeI_4471',snr_min=15,telluric=True)
HD2905=line_study('HD2905','FIES',He_I3,4468,4475,0.015,lim1=1.2,lim2=1.4,plt_mean=True,target='HeI_4471',snr_min=15,telluric=True)

# In[Zero_Moment]
# Zero moment of wind lines: 
folders='4552','H_alpha_prb','H_beta_4861','HeI_5875','4661'
lines=4552.622,6562.80,4861.325,5875.62,4661.6324
bar = progressbar.ProgressBar(max_value=len(lines))
I=0
for i in range(len(lines)):
    try: os.remove('../stars/HD2905/'+folders[i]+'/'+'00_zero_moment.txt')
    except: pass
    SM=reading_summary('../stars/HD2905/'+folders[i]+'/00_HD2905_'+folders[i]+'.txt')
    files=os.listdir('../stars/HD2905/'+folders[i])
    f=open('../stars/HD2905/'+folders[i]+'/'+'00_zero_moment.txt', 'a')
    for fil in range(len(SM.line_file)):
        wv,v,fl=reading_line_file('../stars/HD2905/'+folders[i]+'/'+SM.line_file[fil][1:-1])
        info=SM.Original_file[fil][2:-1],SM.line_file[fil][1:-1],SM.HJD[fil],zero_moment(wv,fl,lines[i])
        f.writelines(str(info))
        f.writelines('\n')
    I+=1
    bar.update(I)
    f.close()

# In[]

SiIII_1='../stars/HD2905/4552/00_HD2905_4552.txt'  
SiIII_2='../stars/HD2905/5739/00_HD2905_5739.txt'

f1,f2 = open(SiIII_1,'r'),open(SiIII_2,'r')
#dat1,dat2 = [],[]
HJD1,fm1,error_fm1=[],[],[]
HJD2,fm2,error_fm2=[],[],[]

HJD1_m,fm1_m,error_fm1_m=[],[],[]
HJD2_m,fm2_m,error_fm2_m=[],[],[]
HJD1_f,fm1_f,error_fm1_f=[],[],[]
HJD2_f,fm2_f,error_fm2_f=[],[],[]
'''
for lin in f1: 
    l = lin.split()  
    try:
        HJD1.append(float(l[1].replace(',',''))),fm1.append(float(l[9].replace(',','')))
        error_fm1.append(float(l[10].replace(',','')))
    except ValueError:
        continue
    
for lin in f2: 
    l = lin.split()
    try:  
        HJD2.append(float(l[1].replace(',',''))),fm2.append(float(l[9].replace(',','')))
        error_fm2.append(float(l[10].replace(',','')))
    except ValueError:
        continue
'''

for lin in f1: 
    l = lin.split() 
    if l[2].replace(',','')=="'MERCATOR'":
        HJD1_m.append(float(l[1].replace(',',''))),fm1_m.append(float(l[9].replace(',','')))
        error_fm1_m.append(float(l[10].replace(',','')))
    elif l[2].replace(',','')=="'FIES'":
        HJD1_f.append(float(l[1].replace(',',''))),fm1_f.append(float(l[9].replace(',','')))
        error_fm1_f.append(float(l[10].replace(',','')))

for lin in f2: 
    l = lin.split() 
    if l[2].replace(',','')=="'MERCATOR'":
        HJD2_m.append(float(l[1].replace(',',''))),fm2_m.append(float(l[9].replace(',','')))
        error_fm2_m.append(float(l[10].replace(',','')))
    elif l[2].replace(',','')=="'FIES'":
        HJD2_f.append(float(l[1].replace(',',''))),fm2_f.append(float(l[9].replace(',','')))
        error_fm2_f.append(float(l[10].replace(',','')))


lis=sorted(zip(HJD1_m,fm1_m,error_fm1_m))
HJD1_m,fm1_m,error_fm1_m=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]
lis=sorted(zip(HJD2_m,fm2_m,error_fm2_m))
HJD2_m,fm2_m,error_fm2_m=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]

lis=sorted(zip(HJD1_f,fm1_f,error_fm1_f))
HJD1_f,fm1_f,error_fm1_f=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]
lis=sorted(zip(HJD2_f,fm2_f,error_fm2_f))
HJD2_f,fm2_f,error_fm2_f=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]



plt.close('all')
plt.errorbar(fm1_m,fm2_m,xerr=error_fm1_m,yerr=error_fm2_m,fmt='.',label='MERCATOR')
plt.errorbar(fm1_f,fm2_f,xerr=error_fm1_f,yerr=error_fm2_f,fmt='.',label='FIES')
plt.plot(fm1,fm2,'.')
plt.grid()
plt.xlabel('Si III (4552)'),plt.ylabel('Si III (5739)')
plt.legend()
plt.savefig('../stars/HD2905/4552_5739.png')


'''
lis=sorted(zip(HJD1,fm1,error_fm1))
HJD1,fm1,error_fm1=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]
lis=sorted(zip(HJD2,fm2,error_fm2))
HJD2,fm2,error_fm2=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]

plt.plot(fm1,fm2,'.')
plt.grid()
plt.xlabel('Si III (4552)'),plt.xlabel('Si III (5739)')

for i in range(0,len(dat)):
    BJD_Sergio.append(float(dat[i][0]))
    fm1_Sergio.append(float(dat[i][2]))
    
fin=open('../stars/HD2905/5739/00_HD2905_5739prueba.txt',"r")
fout = open('../stars/HD2905/5739/00_cleaned.txt',"w+")

for line in fin:
    dat=line.split()
    old=''
    for d in dat:
        if d.startswith("'5739/"):
            pass
        else:
            if old=='':
                old+=d
            else:
                old=old+' '+d
    fout.writelines(str(old))
    fout.writelines('\n')      
fin.close()
fout.close
'''

#%%
H_beta=4861.325 
#prueba_linea('HD2905','FIES',H_beta,4850,4870)
#prueba_linea('HD2905','MERCATOR',H_beta,4850,4870)
HD2905=line_study('HD2905','FIES',H_beta,4850,4870,0.015,lim1=1.35,lim2=1.45,target='H_beta_prb',plt_mean=True)
HD2905=line_study('HD2905','MERCATOR',H_beta,4850,4870,0.015,lim1=1.2,lim2=1.4,target='H_beta_prb',plt_mean=True)#,lim1=0.85,lim2=1,plt_mean=True)
HD2905=line_study('HD2905','SONG',H_beta,4850,4870,0.015,lim1=1.2,lim2=1.4,target='H_beta_prb',plt_mean=True,snr_min=10)

'''
H_beta=4861.32
runcell(0, '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/He_I.py')
HD2905=ls.line_study('HD2905','FIES',H_beta,4850,4870,0.015,lim1=1.35,lim2=1.45,target='H_beta_4861',plt_mean=True)'''
#%%

H_alpha=6562.80
#prueba_linea('HD2905','FIES',H_alpha,6550,6580)
#prueba_linea('HD2905','MERCATOR',H_alpha,6550,6580)
#HD2905=line_study('HD2905','FIES',HeI,5865,5885,0.015,lim1=0.85,lim2=1,plt_mean=True,target='HeI_5875')    
#HD2905=line_study('HD2905','MERCATOR',HeI,5865,5885,0.015,lim1=0.8,lim2=1,plt_mean=True,target='HeI_5875')
HD2905=line_study('HD2905','SONG',H_alpha,6553,6574,0.015,target='H_alpha_prb',plt_mean=True,telluric=True)
HD2905=line_study('HD2905','FIES',H_alpha,6545,6580,0.015,lim1=2,lim2=2.5,target='H_alpha_prb',plt_mean=True,telluric=True)#,lim1=0.85,lim2=1,plt_mean=True)    
HD2905=line_study('HD2905','MERCATOR',H_alpha,6550,6574,0.015,lim1=1.5,lim2=2.,target='H_alpha_prb',plt_mean=True,telluric=True)#,lim1=0.85,lim2=1,plt_mean=True)














