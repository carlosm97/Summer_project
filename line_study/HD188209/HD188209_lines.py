#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:08:25 2021

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

star='HD188209' # O9.5Iab
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
#axs.set_xlim(5887.5,5897.5)
#axs.set_xlim(4000,7000)
axs.set_ylim(0.,1.2)
axs.grid()
plt.show()
# In[Hydrogen]

H_alpha, H_beta=6562.80,4861.325
plt.close('all')
#prueba_linea(star,'FIES',H_alpha,6552,6573,x='vel')
#prueba_linea(star,'MERCATOR',H_alpha,6552,6573,x='vel')
#prueba_linea(star,'SONG',H_alpha,6552,6573,x='vel')

HD188209=line_study(star,'FIES',H_alpha,6552,6573,0.015,lim1=1.1,sig=15,sel='vel',plt_mean=True,telluric=True,target='H_alpha')
HD188209=line_study(star,'MERCATOR',H_alpha,6552,6573,0.015,lim1=1.6,lim2=1.7,plt_mean=True,telluric=True,target='H_alpha')
HD188209=line_study(star,'SONG',H_alpha,6552,6573,0.015,snr_min=15,lim1=1.4,lim2=1.6,plt_mean=True,telluric=True,target='H_alpha')



prueba_linea(star,'FIES',H_beta,4852.5,4871,x='vel')
prueba_linea(star,'MERCATOR',H_beta,4852.5,4871,x='vel')
prueba_linea(star,'SONG',H_beta,4852.5,4871,x='vel')


HD188209=line_study(star,'FIES',H_beta,4852.5,4871,0.015,plt_mean=True,target='H_beta')
HD188209=line_study(star,'MERCATOR',H_beta,4852.5,4871,0.015,telluric=True,plt_mean=True,target='H_beta')
HD188209=line_study(star,'SONG',H_beta,4852.5,4871,0.015,lim1=1.4,lim2=1.5,telluric=True,snr_min=15,plt_mean=True,target='H_beta')
# In[Metals]

plt.close('all')
SiIV=4116.103 #NO SONG
SiIII_1=4552.622
SiIII_2=4567.840
CIV_1=4685.82 #!
OIII_2=4921.96 #!
OIII_1=5592.252 #! 5589.9

#SiIII_3=4574.757 # Con emision, y teniendo las otras... 
#OII_1=4590.974 # Muy al límite... 
#OII_2=4661.6324 # Muy al límite... 
#SiIII_4=5739.734
#SiII_1=6347.11
#CIV_2=5801.35 #!ñ Se mezcla con otra...
#metal_lines=(SiIV,SiIII_1,SiIII_2,OIII_1,OIII_2,CIV_1)

lin=SiIV
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.7,lim2=0.8,normalization=False)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True,lim1=0.8,lim2=0.9,normalization=False)
lin=SiIII_1
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=150,sig=5)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True,lim2=1.1)

lin=SiIII_2
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=120,sig=5)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=120,sig=5)
HD188209=line_study(star,'SONG',lin,lin-5,lin+5,0.015,telluric=True,snr_min=20,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=120,sig=5)

lin=CIV_1
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=0.8,lim2=1,plt_mean=True,sel='vel',sig=10)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True)
HD188209=line_study(star,'SONG',lin,lin-5,lin+5,0.015,lim1=1.15,snr_min=15,telluric=True,plt_mean=True)

lin=OIII_2
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim2=1.1,plt_mean=True)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=1.1,lim2=1.2,telluric=True,plt_mean=True)
HD188209=line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,lim1=1.2,lim2=1.25,sig=5,sel='vel',telluric=True,plt_mean=True)

lin=OIII_1
HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim2=1.1,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=150)
HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=1.1,lim2=1.2,telluric=True,plt_mean=True,automatic_selection=False,sel='vel',ll1=-150,ll2=150)



# In[HeI]
plt.close('all')
HeI_1=4387.929
HeI_2=4713.145
HeI_3=5015.678
HeI_4=5875.62
He=[HeI_1,HeI_2,HeI_3]
# HeI_1 not in SONG 
'''
for lin in He:
    if lin==HeI_1:
        try:
            #line_study(star,'FIES',lin,4383.5,lin+5,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
            line_study(star,'MERCATOR',lin,4383.5,lin+5,0.015,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
            line_study(star,'SONG',lin,4383.5,lin+5,0.015,snr_min=15,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
        except:
            pass
    else:
        try:
            #line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
            line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
            line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4],save=False)
        except:
            pass
'''
try:
    lin=HeI_1
    
    HD188209=line_study(star,'FIES',lin,4383.5,lin+5,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],automatic_selection=False,ll1=-150,ll2=150,sel='vel')
    HD188209=line_study(star,'MERCATOR',lin,4383.5,lin+5,0.015,telluric=True,lim1=0.9,lim2=1,plt_mean=True,target='HeI_'+str(lin)[:4]) 
except:
    pass
try:
    lin=HeI_2
    
    HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],lim1=0.65,sel='vel',sig=5)
    HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim2=1.1,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4]) 
    HD188209=line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,lim1=1.1,sel='vel',sig=5,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4])
except:
    pass
try:        
    lin=HeI_3
    
    HD188209=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],sel='vel',sig=5)
    HD188209=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4])
    HD188209=line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,lim1=1.2,lim2=1.3,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4])  
except:
    pass
try:
    lin=HeI_4
    HD188209=line_study(star,'FIES',lin,lin-6,lin+6,0.015,plt_mean=True,target='HeI_'+str(lin)[:4],sel='vel',sig=5)
    HD188209=line_study(star,'MERCATOR',lin,lin-6,lin+6,0.015,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4])
    HD188209=line_study(star,'SONG',lin,lin-6,lin+6,0.015,snr_min=15,sel='vel',sig=5,telluric=True,plt_mean=True,target='HeI_'+str(lin)[:4])
except:
    pass
#HD188209=line_study(star,'SONG',lin,lin-4,lin+5,0.015,snr_min=15,telluric=True,plt_mean=True,save=False,target='HeI_'+str(lin)[:4])


# In[HeII]
plt.close('all')
HeII_1=4541.591
HeII_2=4859.319
HeII_3=5411.52
He=[HeII_1,HeII_2,HeII_3]

lin=HeII_1
HD188209=line_study(star,'FIES',lin,lin-6,lin+6,0.015,plt_mean=True,normalization=False,target='HeII_'+str(lin)[:4],\
                    automatic_selection=False,sel='vel',ll1=-200,ll2=150,sig=5)
HD188209=line_study(star,'MERCATOR',lin,lin-6,lin+6,0.015,telluric=True,plt_mean=True,target='HeII_'+str(lin)[:4],\
                    automatic_selection=False,sel='vel',ll1=-200,ll2=150,sig=5,normalization=False)
HD188209=line_study(star,'SONG',lin,lin-6,lin+6,0.015,telluric=True,plt_mean=True,target='HeII_'+str(lin)[:4],\
                    automatic_selection=False,sel='vel',ll1=-200,ll2=150,sig=5,snr_min=15)

lin=HeII_2
HD188209=line_study(star,'FIES',lin,lin-6,lin+7,0.015,plt_mean=True,target='HeII_'+str(lin)[:4])
HD188209=line_study(star,'MERCATOR',lin,lin-6,lin+7,0.015,plt_mean=True,target='HeII_'+str(lin)[:4],telluric=True)
HD188209=line_study(star,'SONG',lin,lin-6,lin+7,0.015,plt_mean=True,target='HeII_'+str(lin)[:4],telluric=True,snr_min=15,sel='vel',sig=5)

lin=HeII_3
HD188209=line_study(star,'FIES',lin,lin-6,lin+6,0.015,lim1=0.8,sel='vel',sig=5,plt_mean=True,target='HeII_'+str(lin)[:4])
HD188209=line_study(star,'MERCATOR',lin,lin-6,lin+6,0.015,lim1=0.95,lim2=1,telluric=True,plt_mean=True,target='HeII_'+str(lin)[:4])
HD188209=line_study(star,'SONG',lin,lin-6,lin+6,0.015,snr_min=15,telluric=True,plt_mean=True,target='HeII_'+str(lin)[:4],\
                    automatic_selection=False,sel='vel',ll1=-150,ll2=190,sig=5)
   





































