#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:26:36 2021

@author: charlie
"""

# Vamos a añadir FIES. 

import numpy as np 
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy import integrate
import formulas as fmr
from scipy.interpolate import interp1d
from db import *
from spec import *
from RV import *

HD2905_spectrum=findstar('HD2905')
try:
    del HD2905_MERCATOR
except NameError:
    print('No previous HD2905_MERCATOR')
try:
    del HD2905_FIES
except NameError:
    print('No previous HD2905_FIES')
i=0
HD2905_MERCATOR=[MERCATOR for MERCATOR in HD2905_spectrum if '_M_' in MERCATOR]
HD2905_FIES=[FIES for FIES in HD2905_spectrum if '_N_' in FIES]


def zero_moment(x_l,y_l,line):
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int0=integrate.cumtrapz(F_l, v_l, initial=0)    
    return y_int0[-1]    

def first_moment(x_l,y_l,line):
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int1 = integrate.cumtrapz(F_l*v_l, v_l, initial=0)
    return y_int1[-1]




# In[]

# Interpolate a function in the same wavelengths all of them. 
try:
    del X_int,Y_int,HJD,RV,snr
except NameError:
    print('Ignored')

SiIII,c=4552.662, 299792.458 
X_int,Y_int,HJD,RV,snr=[],[],[],[],[]

for sMerc in HD2905_MERCATOR:
    rv0=RV0(SiIII,sMerc)
    SM=spec(sMerc)
    #wf=SM.waveflux(4549, 4556)
    SM.resamp(0.015,lwl=4549,rwl=4556)
    X_int.append(SM.wave),Y_int.append(SM.flux),HJD.append(SM.hjd),RV.append(RV0(SiIII,sMerc))
    #SM=spec(sMerc)
    #SM.resamp(0.015,lwl=4549,rwl=4556) 
    #X_int.append(SM.wave),Y_int.append(SM.flux),HJD.append(SM.hjd),RV.append(RV0(4552.622,sMerc))
    snr.append(SM.snr)

# In[]
# Computing the mean spectrum. 
plt.close('all')
X_int,Y_int=np.vstack(X_int),np.vstack(Y_int)
X_mean,Y_mean=np.mean(X_int,0),np.mean(Y_int,0)
Error=np.nanstd(Y_int,0)
v_x=(X_mean-SiIII)*c/SiIII 
plt.errorbar(v_x,Y_mean,yerr=Error,fmt='.',color='k',label='Mean spectrum')
[plt.plot(v_x,Y_int[i,:],alpha=0.3) for i in range(0,len(X_int))]
plt.xlabel('Velocity [km/s]'),plt.ylabel('Normalized flux'),plt.title('SiIII')


plt.grid()
plt.ylim(0.7,1.05),plt.xlim(v_x[0],v_x[-1])
plt.legend()

# In[]
# Computing the TSV and deciding sigma. 

std1=[np.std(Y_int[i][:60]) for i in range(0,len(X_int))]
std2=[np.std(Y_int[i][-60:]) for i in range(0,len(X_int))]
std=(np.mean(std1)+np.mean(std2))/2

plt.figure()
plt.plot(v_x,Error)
plt.axhline(y=1.2*np.mean(std),color='g',label='Cut criteria (1.2 mean std)')
plt.axhline(y=1.*np.mean(std),color='r',label='Cut criteria (mean std)')
plt.ylabel('TVS')
plt.xlabel('Velocity [km/s]')
plt.legend()
plt.grid()
xg,yg=[],[]
M = int(len(Error)/2)


for i in range(0,M):
    mmax_1p1=X_mean[M+i]
    if Error[M+i]>1.2*np.mean(std):
        continue
    else: 
        break
for i in range(0,M):
    mmin_1p1=X_mean[M-i]
    if Error[M-i]>1.2*np.mean(std):
        continue
    else: 
        break
    
for i in range(0,M):
    mmax_1=X_mean[M+i]
    if Error[M+i]>1.*np.mean(std):
        continue
    else: 
        break
for i in range(0,M):
    mmin_1=X_mean[M-i]
    if Error[M-i]>1.*np.mean(std):
        continue
    else: 
        break    
# Vamos a escoger un criterio para sigma tal que la FWHM de la distribución sea 
# justamente la correspondiente a 1. 
sigma=((mmin_1p1-mmin_1)+(mmax_1-mmax_1p1))/2/2.355
# Pasarlo a velocidades para plotear

plt.axvline(x=(mmin_1-sigma*2.355-SiIII)*c/SiIII,color='g',label='Cut criteria (1.1 mean std)')
plt.axvline(x=(mmin_1+sigma*2.355-SiIII)*c/SiIII,color='g',label='Cut criteria (1.1 mean std)')

plt.axvline(x=(mmax_1-sigma*2.355-SiIII)*c/SiIII,color='g',label='Cut criteria (1.1 mean std)')
plt.axvline(x=(mmax_1+sigma*2.355-SiIII)*c/SiIII,color='g',label='Cut criteria (1.1 mean std)')


# In[]
# La RV sacada por el método del programa, tiene un sesgo mayor, así que ajustamos
# una gaussiana. 
# Buscamos sacar los errores por un MC
def fit_gaussian(central_ojo,ristra,pars0=None):
    x,y=ristra[0],ristra[1]
    #mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))#y[E_min[1]-400:E_min[1]+400]<=0.93    
    # Ajuste de gaussiana al perfil 
    X_E=[x,y]#[x_E[x_E !=0],y_E[y_E !=0]]
    if pars0==None:
        pars0=Pars0(X_E)
    
    pairs, covar=fmr.adjust(fmr.gaussian_flux,X_E,pars0)
    x_lin=x#np.linspace(x_E[x_E !=0][0],x_E[x_E !=0][-1],100)
    return [x_lin,fmr.gaussian_flux(x_lin,*pairs),pairs]
def radial_velocity(x_l,y_l,SiIII,pars0=None):
    v_l = (x_l-SiIII)*c/SiIII    
    return fit_gaussian(0,[v_l,y_l],pars0)[-1][-1]


first_m,error_fm=[],[]
rv,error_rv=[],[]
I=0
BJD_list=[]

for j in range(0,len(X_int)):
    v,mini=[],[]
    #BJD,std,x_Si,y_Si=fmr.get_data('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/SiIII_prueba/'+j)
    #BJD_list.append(BJD)
    #Z=sorted(zip(BJD_list,m1_def,mini))
    #BJD_list,m1_def,mini=[Z[i][0] for i in range(0,len(BJD_list))],[Z[i][1] for i in range(0,len(BJD_list))],[Z[i][2] for i in range(0,len(BJD_list))]
    x_Si,y_Si= X_int[j], Y_int[j]
    for i in range(0,500):
        rmmin=np.random.normal(loc=mmin_1,scale=2*sigma)
        rmmax=np.random.normal(loc=mmax_1,scale=2*sigma)
        ran1,ran2=fmr.find_nearest(x_Si,rmmin)[1],fmr.find_nearest(x_Si,rmmax)[1]
        m0=zero_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],SiIII)
        m1=first_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],SiIII)
        v.append(m1/m0)
       
        mini.append(radial_velocity(x_Si[ran1:ran2],y_Si[ran1:ran2],SiIII,[0.2,5,0]))
    #plt.figure()     
    n, bins = np.histogram(v, bins=100,density='true')
    mids = 0.5*(bins[1:] + bins[:-1])
    mean = np.average(mids, weights=n)
    var = np.average((mids - mean)**2, weights=n)
    error = np.sqrt(var)
    first_m.append(mean),error_fm.append(error)   
    
    n, bins = np.histogram(mini, bins=100,density='true')
    mids = 0.5*(bins[1:] + bins[:-1])
    mean = np.average(mids, weights=n)
    var = np.average((mids - mean)**2, weights=n)
    error = np.sqrt(var)
    rv.append(mean),error_rv.append(error)
    I=I+1
    print(I)
bar.finish()
plt.figure()
plt.errorbar(np.array(HJD)-2450000,first_m,yerr=error_fm,fmt='.',label='First moment')
plt.errorbar(np.array(HJD)-2450000,rv,yerr=error_rv,fmt='.',label='Gaussian radial velocity')
plt.plot(np.array(HJD)-2450000,RV,'.',label='Program radial velocity')
plt.grid()
plt.tick_params(axis='x', labelsize=17),plt.tick_params(axis='y', labelsize=17)
plt.ylabel('Velocity',fontsize=17), plt.xlabel(r'BJD-2450000 [days]',fontsize=17)
plt.legend()
# In[]
plt.close('all')
plt.figure()
desv_rv=np.std(np.array(rv)-np.array(RV)) # Generamos un error tal que ambas medidas
# de la velocidad radial sean consistentes. 
# El error así calculado es mayor al calculado a través de MC. Siendo conservadores, 
# parece mejor decisión tomar éste. 
plt.errorbar(RV,rv,yerr=desv_rv,xerr=desv_rv,fmt='.'),plt.grid()

plt.plot([min(RV),max(RV)],[min(RV),max(RV)])
plt.xlabel('Program radial velocity'),plt.ylabel('Gaussian radial velocity')

# Las velocidades calculadas de ambas formas es muy consistente: forma una recta 
# muy clara. Por un lado, al calcularla como gaussiana, permite hacer MC y estimar 
# el error. El sacarlo con pyIACOB, sin embargo, es mucho más rápido, y el error 
# parece lo suficientemente despreciable como para que no merezca la pena sacarlo.


# In[]
plt.close('all')
plt.figure()
#plt.errorbar(first_m_g,rv_g,xerr=error_fm_g,yerr=error_rv_g,fmt='.',label='Gaussian')
plt.errorbar(first_m,rv,xerr=error_fm,yerr=desv_rv,fmt='.',label='TVS')
plt.xlabel('First moment [km/s]')
plt.ylabel('Radial velocity [km/s]')
plt.plot([min(rv),max(rv)],[min(rv),max(rv)],label='1-1 line')
plt.grid()   

linear_model=np.polyfit(first_m,rv,1)
linear_model_fn=np.poly1d(linear_model)
pendiente=linear_model[0]
D_pendiente=np.sqrt(len(first_m))*max(error_fm)/np.sqrt(len(first_m)*sum(np.array(first_m)**2)-(sum(first_m))**2)
print('Slope=',pendiente,r'$\pm$',D_pendiente)
plt.plot(first_m,linear_model_fn(first_m),color="green", label='Fit line')
plt.title('Slope='+str(round(pendiente,3))+r'$\pm$'+str(round(D_pendiente,3)))
plt.legend()
plt.show()

# In[]

# Guardamos las líneas en ficheros: 
import csv
    
for i in range(0,len(HD2905_MERCATOR)):
    info=zip(X_int[i],(X_int[i]-SiIII)*c/SiIII,Y_int[i])
    BJD_s=str(HJD[i])
    BJD_s=BJD_s.replace('.','p')
    with open('../stars/HD2905/4552/HD_2905_'+BJD_s+'_MERCATOR_4552.txt', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(info)
    f.close()
# nombre, JD, intrumento, resolución, mínimo y máximo de la línea, sigma a su alrededor, 
# 
f=open('../stars/HD2905/4552/00_HD_2905_.txt', 'w')
f.writelines('File,     HJD,    Instrument,     Resolution,     lam_min,    lam_max,\
    sigma,     snr,    first moment,   error fm,   Radial velocity,    error rv, \
    Other error rv')
f.writelines('\n')                
for i in range(0,len(HD2905_MERCATOR)):
    info=HD2905_MERCATOR[i].replace('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/HD2905_FH/',''),\
        HJD[i],'MERCATOR',85000,mmin_1,mmax_1,sigma,snr[i],first_m[i],error_fm[i],\
        rv[i],error_rv[i],desv_rv
    f.writelines(str(info))
    f.writelines('\n')
f.close()            

# In[TVS_example]
plt.close('all')
spectrum=findstar(star)
line,step=4552.622,0.015
ll,rl=line-6,line+6
star='HD2905'
c1,c2=50,-50
condition='_N_'  
ext='_FH'  
files=[ins for ins in spectrum if condition in ins]
snr_min=0
lim1=1
automatic_selection=False
# Interpolate a function in the same wavelengths all of them. 
try:
    del X_int,Y_int,HJD,RV
except NameError:
    pass
c=299792.458 
X_int,Y_int,HJD,RV=[],[],[],[]
bar = progressbar.ProgressBar(max_value=len(files))
I=0
ll,rl=line-step*int((line-ll)/step),line+step*int((rl-line)/step)
telluric=False
normalization=True
fig,axs=plt.subplots(2,sharex=True,gridspec_kw={'height_ratios': [1,2]})
for s in files:
    #rv0=RV0(line,s,verbose=False,line=line)
    try:
        SM=sp.spec(s,line=line)
        rv0=SM.fitline(line)['RV_kms']
        SM.resamp(step,lwl=ll,rwl=rl)
        wv,fl=SM.wave,SM.flux
        if normalization==True:
            pol=np.polyfit(np.concatenate((wv[:c1],wv[c2:])),np.concatenate((fl[:c1],fl[c2:])),1)
            corr=fl/np.polyval(pol,wv)   
        else:
            corr=fl
        if 1/np.nanstd(np.concatenate((corr[:c1],corr[c2:])))>snr_min: # We avoid data with really bad SNR
            if telluric==True:
                CRR=cosmic_telluric(corr)
            else: 
                CRR=corr
            X_int.append(wv),Y_int.append(CRR),HJD.append(SM.hjd),RV.append(rv0)
        v_x=(wv-line)*c/line
        if np.random.rand()>0.5:
            axs[1].plot(v_x,Y_int[-1])
    except: 
        pass
    I=I+1
    bar.update(I)
#[print(s) for s in files]
# Computing the mean spectrum. 
X_int,Y_int=np.vstack(X_int),np.vstack(Y_int)
X_mean,Y_mean=np.mean(X_int,0),np.mean(Y_int,0)
Error=np.nanstd(Y_int,0)
v_x=(X_mean-line)*c/line 
axs[1].plot(v_x,Y_mean,'k')
# Computing the TSV and deciding sigma. 
std1=[np.std(Y_int[i][:c1]) for i in range(0,len(X_int))]
std2=[np.std(Y_int[i][c2:]) for i in range(0,len(X_int))]
STD=(np.array(std1)+np.array(std2))/2
std=(np.mean(std1)+np.mean(std2))/2
axs[0].plot(v_x,Error)
axs[0].axhline(y=lim1*np.mean(std),color='g',label='Cut criteria ('+str(lim1)+' mean std)')
axs[0].set_ylabel('TVS')
axs[1].set_xlabel('Velocity [km/s]'),axs[1].set_ylabel('Normalized flux')
#axs[0].legend()
axs[0].grid()
xg,yg=[],[]
M = int(len(Error)/2)
for i in range(0,M):
    mmax_1=X_mean[M+i]
    if Error[M+i]>lim1*np.mean(std):
        continue
    else: 
        break
for i in range(0,M):
    mmin_1=X_mean[M-i]
    if Error[M-i]>lim1*np.mean(std):
        continue
    else: 
        break    

#axs[0].title(str(line))
axs[1].plot(v_x,Y_mean)
axs[1].grid()
x=v_x[Error<np.mean(std)]
max(x[x<0])
axs[0].axvline(x=min(x[x>0]),color='r',label='Cut criteria (1.1 mean std)')
axs[0].axvline(x=max(x[x<0]),color='r',label='Cut criteria (1.1 mean std)')
axs[1].axvline(x=min(x[x>0]),color='r',label='Cut criteria (1.1 mean std)')
axs[1].axvline(x=max(x[x<0]),color='r',label='Cut criteria (1.1 mean std)')


#    plt.savefig('../stars/'+star+'/'+target+'/TVS_criteria_'+instrument+'.png')
















