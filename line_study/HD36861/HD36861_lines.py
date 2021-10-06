#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 14:49:40 2021

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

star='HD36861' # O8III
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

prueba_linea(star,'FIES',H_alpha,6540,6585)

HD91316=line_study(star,'FIES',H_beta,4845,4880,0.015,plt_mean=True,target='H_beta',sig=30,automatic_selection=False,ll1=-700,ll2=700,sel='vel')
HD91316=line_study(star,'MERCATOR',H_beta,4845,4880,0.015,plt_mean=True,target='H_beta',sig=30,automatic_selection=False,telluric=True,ll1=-750,ll2=700,sel='vel')
HD91316=line_study(star,'SONG',H_beta,4850,4873,0.015,plt_mean=True,target='H_beta',snr_min=15,telluric=True,sig=30,automatic_selection=False,ll1=-500,ll2=500,sel='vel')

HD91316=line_study(star,'FIES',H_alpha,6550,6575,0.015,plt_mean=True,telluric=True,target='H_alpha',sig=10,automatic_selection=False,ll1=-400,ll2=400,sel='vel')
HD91316=line_study(star,'MERCATOR',H_alpha,6550,6575,0.015,plt_mean=True,telluric=True,target='H_alpha',sig=10,automatic_selection=False,ll1=-400,ll2=400,sel='vel')
HD91316=line_study(star,'SONG',H_alpha,6552,6574,0.015,plt_mean=True,telluric=True,target='H_alpha',snr_min=15,sig=10,automatic_selection=False,ll1=-400,ll2=400,sel='vel')

# In[Metals_try]
SiIV=4116.103
# SiIII_1=4552.622
# SiIII_2=4567.840
# SiIII_3=4574.757
# OII_1=4590.974
# OII_2=4661.6324
# SiIII_4=5739.734
OIII_1=5592.252
# SiII_1=6347.11

prueba_linea(star,'FIES',SiIV,4110,4120)
# prueba_linea(star,'FIES',SiIII_1,4545,4560) # Floja
# prueba_linea(star,'FIES',SiIII_2,4560,4575) # 4568¿?
# prueba_linea(star,'FIES',SiIII_3,4570,4580) # Not
# prueba_linea(star,'FIES',OII_1,4585,4595) # Not
# prueba_linea(star,'FIES',OII_2,4655,4665) # Not
# prueba_linea(star,'FIES',SiIII_4,5735,5745) # Not
prueba_linea(star,'FIES',OIII_1,5585,5600)
# prueba_linea(star,'FIES',SiII_1,6340,6355) # Not

# SII_1=5055.984
# SII_2=5453.855
# SII_3=5639.977
# SiII_2=6347.11
# NeI=6402.2484
# CII_1=4267.183
# CII_2=6582.88

# prueba_linea(star,'FIES',SII_1,5050,5060) # Not
# prueba_linea(star,'FIES',SII_2,5445,5460) # Not
# prueba_linea(star,'FIES',SII_3,5635,5645) # Not
# prueba_linea(star,'FIES',SiII_2,6340,6355) # Not
# prueba_linea(star,'FIES',NeI,6395,6410) # Not
# prueba_linea(star,'FIES',CII_1,4260,4275) # Not
# prueba_linea(star,'FIES',CII_2,6575,6590)  # Not


# In[Metals]
SiIV=4116.103
OIII_1=5592.252 #!
CIV_1=4685.82 #!
OIII_2=4921.96 #!
CIV_2=5801.35 #!
# FeIII=5812.200 # Waiting for Sergio... 

HD36861=line_study(star,'FIES',SiIV,4112,4121,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=170,sel='vel',sig=3)
HD36861=line_study(star,'MERCATOR',SiIV,4112,4121,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=170,sel='vel',sig=3,telluric=True)

HD36861=line_study(star,'FIES',OIII_1,5588,5598,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=200,sel='vel',sig=5)
HD36861=line_study(star,'MERCATOR',OIII_1,5588,5598,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=200,sel='vel',sig=5,telluric=True)
HD36861=line_study(star,'SONG',OIII_1,5590.3,5596,0.015,plt_mean=True,automatic_selection=False,ll1=-75,ll2=175,sel='vel',sig=5,telluric=True)

HD36861=line_study(star,'FIES',CIV_1,4680,4693,0.015,plt_mean=True,automatic_selection=False,ll1=-300,ll2=300,sel='vel',sig=10)
HD36861=line_study(star,'MERCATOR',CIV_1,4680,4693,0.015,plt_mean=True,automatic_selection=False,ll1=-200,ll2=250,sel='vel',sig=10,telluric=True)
HD36861=line_study(star,'SONG',CIV_1,4680,4691.5,0.015,plt_mean=True,automatic_selection=False,ll1=-150,ll2=250,sel='vel',sig=10,telluric=True)


HD36861=line_study(star,'FIES',OIII_2,4917,4927,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'MERCATOR',OIII_2,4917,4927,0.015,telluric=True,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'SONG',OIII_2,4917,4927,0.015,telluric=True,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')


HD36861=line_study(star,'FIES',CIV_2,5796,5806,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'MERCATOR',CIV_2,5796,5806,0.015,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel',telluric=True)
# No in SONG

# In[HeI]
HeI_1=4387.929
HeI_2=4713.145
HeI_3=5015.678
HeI_4=5875.62


HD36861=line_study(star,'FIES',HeI_1,4383.6,4393,0.015,target='HeI_4387',plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'MERCATOR',HeI_1,4383.6,4393,0.015,target='HeI_4387',telluric=True,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'SONG',HeI_1,4383.6,4393,0.015,target='HeI_4387',telluric=True,snr_min=15,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')

HD36861=line_study(star,'FIES',HeI_2,4708,4718,0.015,target='HeI_4713',plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'MERCATOR',HeI_2,4708,4718,0.015,target='HeI_4713',telluric=True,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')
HD36861=line_study(star,'SONG',HeI_2,4708,4718,0.015,target='HeI_4713',telluric=True,snr_min=15,plt_mean=True,automatic_selection=False,ll1=-100,ll2=150,sel='vel')

HD36861=line_study(star,'FIES',HeI_3,5012,5020,0.015,target='HeI_5015',plt_mean=True,automatic_selection=False,ll1=-80,ll2=150,sel='vel')
HD36861=line_study(star,'MERCATOR',HeI_3,5012,5020,0.015,target='HeI_5015',telluric=True,plt_mean=True,automatic_selection=False,ll1=-80,ll2=150,sel='vel')
HD36861=line_study(star,'SONG',HeI_3,5012,5020,0.015,target='HeI_5015',snr_min=15,plt_mean=True,automatic_selection=False,ll1=-80,ll2=150,sel='vel')

HD36861=line_study(star,'FIES',HeI_4,5872,5881,0.015,target='HeI_5875',plt_mean=True,automatic_selection=False,ll1=-100,ll2=180,sel='vel')
HD36861=line_study(star,'MERCATOR',HeI_4,5872,5881,0.015,target='HeI_5875',telluric=True,plt_mean=True,automatic_selection=False,ll1=-100,ll2=180,sel='vel')
HD36861=line_study(star,'SONG',HeI_4,5872,5881,0.015,target='HeI_5875',plt_mean=True,automatic_selection=False,ll1=-100,ll2=180,sel='vel')

# In[HeII]
HeII_1=4541.591
HeII_2=4859.319
HeII_3=5411.52


HD36861=line_study(star,'FIES',HeII_1,4537,4547,0.015,target='HeII_4541',plt_mean=True,automatic_selection=False,ll1=-200,ll2=200,sel='vel',sig=5)
HD36861=line_study(star,'MERCATOR',HeII_1,4537,4547,0.015,target='HeII_4541',telluric=True,plt_mean=True,automatic_selection=False,ll1=-200,ll2=200,sel='vel',sig=5)
HD36861=line_study(star,'SONG',HeII_1,4537,4547,0.015,target='HeII_4541',telluric=True,snr_min=15,plt_mean=True,automatic_selection=False,ll1=-200,ll2=200,sel='vel',sig=5)


HD36861=line_study(star,'FIES',HeII_2,4853,4869,0.015,target='HeII_4859',plt_mean=True,automatic_selection=False,ll1=-300,ll2=400,sel='vel',sig=8)
HD36861=line_study(star,'MERCATOR',HeII_2,4853,4869,0.015,target='HeII_4859',telluric=True,plt_mean=True,automatic_selection=False,ll1=-300,ll2=400,sel='vel',sig=8)
HD36861=line_study(star,'SONG',HeII_2,4853,4869,0.015,target='HeII_4859',plt_mean=True,automatic_selection=False,ll1=-300,ll2=400,sel='vel',sig=8)

