#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:35:13 2021

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
import spec as sp
from RV import *
import line_study as ls
import formulas as fmr
import progressbar
#%%
path='../stars/HD2905/H_alpha_6562/'
dre=os.listdir(path[:-1])
fil_flux=[]
for file in dre:
    if file.endswith('.txt') and not file.startswith('00') and 'SONG' in file:
        wv,v,f=ls.reading_line_file(file,path=path)
        fil_flux.append(f)
#plt.plot(wv,f)
plt.imshow(fil_flux)

#%%
path='../stars/HD2905/H_beta_4861/'
dre=os.listdir(path[:-1])
fil_flux=[]
for file in dre:
    if file.endswith('.txt') and not file.startswith('00') and 'MERCATOR' in file:
        wv,v,f=ls.reading_line_file(file,path=path)
        fil_flux.append(f)
#plt.plot(wv,f)
plt.imshow(fil_flux)

#%%
'''
iter=3
t=3
plt.close('all')
bar = progressbar.ProgressBar(max_value=t*len(fil_flux))
I=0
for rep in range(t):
    D=fil_flux
    # t== Último autovector que considero. 
    # Saco la matriz de correlación
    C=np.matrix(D).transpose()*np.matrix(D)
    #Diagonalizamos y ordenamos 
    evals, evecs = np.linalg.eig(C)
    B=evecs.copy()
    b=np.array(B)
    # Sacamos la matriz de coeficientes en la nueva base, diagonalizada 
    A=np.dot(D,b)#[:,:,0])
    D_PCA=np.real(np.dot(A[:,:t],np.matrix(B[:,:t]).H))   
    for f in range(len(fil_flux)):
        #fig,(ax1,ax2,ax3)=plt.subplots()#, sharex=True)
        #ax1.imshow(D)
        #ax2.imshow(D_PCA)
        #ax3.plot(wv,fil_flux[3],label='Original')
        #ax3.plot(wv,D_PCA[3].H,label='Reduced') 
        
        #plt.figure()
        #plt.plot(wv,(fil_flux/D_PCA)[3].H,'x')
        std=np.std((fil_flux/D_PCA)[f].H)
        #plt.axhline(y=1-2*std)
        #plt.grid()
        # Con esto tenemos localizados posibles candidatos a líneas telúricas. Éstas 
        # deberías encajas bien con gaussianas, así que podemos "modelarlas". 
        save=[]
        j=-1
        I+=1
        bar.update(I)   
        for i in range(1,np.shape(fil_flux)[1]):
            if (fil_flux/D_PCA)[f,i]<1-2*std and (fil_flux/D_PCA)[f,i-1]<1-2*std:
                if (fil_flux/D_PCA)[f,i-2]>=1-2*std:
                    j+=1
                    save.append([])
                    save[j].append((wv[i-1],(fil_flux/D_PCA)[f,i-1]))
                    save[j].append((wv[i],(fil_flux/D_PCA)[f,i]))
                else:
                    save[j].append((wv[i],(fil_flux/D_PCA)[f,i]))
        #save=np.vstack(save)
        wv=np.around(wv,3)
        del_wave=[]
        for j in range(len(save)):
            r=save[j]
            r=np.around(r,3)
            r=np.vstack(r)
            wve=np.array((r[:,0],r[:,1]))
            l1,l2=int(np.where(wv==r[0][0])[0]),int(np.where(wv==r[-1][0])[0])
            #print('initial',l1,l2)
            #l1-=1
            wve,flx=list(r[:,0]),list(r[:,1])
            if l2-l1>=3:
                #wve.insert(0,wv[l1]),flx.insert(0,(fil_flux/D_PCA)[3,l1])
                #while sum((ls.fit_gaussian(r[:,0][r[:,1].argmin()],np.array((wve,flx)))[1]-flx)**2)/len(flx)<=\
                #    1.5*sum((ls.fit_gaussian(r[:,0][r[:,1].argmin()],np.array((wve[1:],flx[1:])))[1]-flx[1:])**2)/len(flx[1:]):
                #        l1-=1
                #        wve.insert(0,wv[l1]),flx.insert(0,(fil_flux/D_PCA)[3,l1])
                #while sum((ls.fit_gaussian(r[:,0][r[:,1].argmin()],np.array((wve,flx)))[1]-flx)**2)/len(flx)<=\
                #    1.5*sum((ls.fit_gaussian(r[:,0][r[:,1].argmin()],np.array((wve[:-1],flx[:-1])))[1]-flx[:-1])**2)/len(flx[:-1]):
                #        l2+=1
                #        wve.append(wv[l2]),flx.append((fil_flux/D_PCA)[3,l2])
                #print('final',l1,l2)    
                while np.mean((fil_flux/D_PCA)[f,l1-2:l1])-np.mean((fil_flux/D_PCA)[f,l1:l1+2])>0:
                    l1-=1
                while np.mean((fil_flux/D_PCA)[f,l2:l2+2])-np.mean((fil_flux/D_PCA)[f,l2-2:l2])>0:
                    l2+=1       
                if wv[l2]-wv[l1]<1:
                    flx=[(fil_flux/D_PCA)[f,k] for k in range(l1,l2)]
                    #ft=ls.fit_gaussian(r[:,0][r[:,1].argmin()],np.array((wv[l1:l2],flx)))
                    plt.plot(wv[l1:l2],fil_flux[f][l1:l2])
                    del_wave.append(wv[l1:l2])
        try:
            del_wave=np.concatenate(del_wave)
        except:
            continue
        wv_no,flx_no=wv[~np.isin(wv,del_wave)],fil_flux[f][~np.isin(wv,del_wave)]
        inter=interp1d(wv_no,flx_no)
        flx_n=inter(wv)
        #plt.plot(wv,flx_n)
        fil_flux[f]=flx_n
        #for j in range(len(fil_flux)):
        #    fil_flux[j]=fil_flux[j][~np.isin(wv,wve)]

    
    

y=(fil_flux/D_PCA)[3].H
wv,sv=np.around(wv,3),np.around(np.vstack(save)[:,0],3)
wv[~np.isin(wv,sv)],fil_flux[3][~np.isin(wv,sv)]

#%%
plt.close('all')
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq, ifft
n = len(fil_flux[3])
Y = fft(fil_flux[3])/n
dl = np.mean(np.diff(wv))
frq = fftfreq(n, dl)

fg,ax=plt.subplots()
ax.plot(frq,Y,'.')
ax.grid()

plt.figure()
Y_int=Y[abs(Y)>np.mean(Y)]
y=ifft(Y_int)
#df = np.mean(np.diff(frq[Y>np.mean(Y)]))
#f = fftfreq(n, df)
plt.plot(y)


plt.figure()
plt.plot(wv,fil_flux[3])
'''
#%%

def cosmic_telluric(flx, sigclip=1.5, iter=3):
    '''
    Function to remove cosmic rays and telluric lines in the spectra.

    Parameters
    ----------
    
    flx : array or array_like, 
        Array containing the original flux withour cleaning. 
    
    sigclip : float, optional
        Sigma clipping value used to remove rays. Default is 1.5.

    iter : int, optional
        Number of iterations of the sigma clipping to remove cosmic rays.

    Returns
    -------
    Corrected flux, without 
    '''
    #www.towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22
    # First we calculated ∇x(i):
    delta_flux = [flx[i+1] - flx[i] for i in np.arange(len(flx) - 1)]

    median_int = np.median(delta_flux)
    mad_int = np.median([np.abs(delta_flux - median_int)])
    modified_z_scores = 0.6745 * (delta_flux - median_int) / mad_int
    # The multiplier 0.6745 is the 0.75th quartile of the standard normal
    # distribution, to which the median absolute deviation converges to.
    flux_norm =  np.concatenate(([1], abs(modified_z_scores))) 
    for i in range(iter):
        std = np.nanstd(flux_norm)
        flux_norm = np.where(abs(flux_norm - 1) > sigclip*std, np.nan, flux_norm) 
    flux_clean = np.where(np.isnan(flux_norm), np.nan, flx)
    #flx_nn,wv_nn=flux_clean[~np.isnan(flux_clean)],wv[~np.isnan(flux_clean)]
    n_p=30
    for j in range(n_p,len(flux_clean) - n_p):
        a=(np.nanmean(flux_clean[j+1:j+n_p])-np.nanmean(flux_clean[j-n_p+1:j]))/(2*n_p-2)
        b=(-a*(j-n_p/2)+np.nanmean(flux_clean[j-n_p+1:j])-a*(j+n_p/2)+np.nanmean(flux_clean[j+1:j+n_p]))/2                                            
        if abs(flux_clean[j]-(a*j+b))>2*(np.nanstd(flux_clean[j-n_p+1:j]/(a*np.linspace(j-n_p+1,j-1,n_p-1)+b))\
         +np.nanstd(flux_clean[j+1:j+n_p]/(a*np.linspace(j+1,j+n_p-1,n_p-1)+b)))/2:
            flux_clean[j]=np.nan
            
    nans = np.isnan(flux_clean)
    x = lambda z: z.nonzero()[0]
    flux_clean[nans] = np.interp(x(nans), x(~nans), flux_clean[~nans])
    flux_clean = np.where((flux_clean > flx+0.01) | (abs(flux_clean-flx) > 0.05),
        flux_clean,flx)
    flux=flux_clean
    return flux
cl_flx=[]
bar = progressbar.ProgressBar(max_value=len(fil_flux))
I=0
for l in range(len(fil_flux)):
    cleaned=cosmic_telluric(fil_flux[l])
    cl_flx.append(cleaned)
    I+=1
    bar.update(I)
#%%
plt.close('all')
fig,axs=plt.subplots(2,2)
'''
corr=fil_flux[190]
axs[0].plot(corr)
axs[1].plot(telluric(corr))
'''
im1=axs[0,0].imshow(fil_flux,vmin=0.8,vmax=1.5)
axs[1,0].imshow(cl_flx,vmin=0.8,vmax=1.5)
cbar_ax = fig.add_axes([0.44, 0.15, 0.05, 0.5])
fig.colorbar(im1, cax=cbar_ax)
axs[1,1].set_xlabel(r'$\lambda [\AA]$')
axs[0,0].set_title('Original'),axs[1,0].set_title('Corrected')
#axs[0].grid(),axs[1].grid()

plt.show()
fi,ax=plt.subplots()
for l in range(len(fil_flux)):
    axs[0,1].plot(wv,fil_flux[l],alpha=0.3)
    axs[1,1].plot(wv,cl_flx[l],alpha=0.3) 
    ax.plot(wv,SM.cosmic)
axs[0,1].plot(wv,np.mean(fil_flux,0),'k')
axs[1,1].plot(wv,np.mean(cl_flx,0),'k')
#%%

spectrum=findstar('HD2905')

condition='_M_'
ext='_FH'
files=[ins for ins in spectrum if condition in ins]
# Interpolate a function in the same wavelengths all of them. 
try:
    del X_int,Y_int,HJD,RV
except NameError:
    pass
c=299792.458 
X_int,Y_int,HJD,RV=[],[],[],[]
bar = progressbar.ProgressBar(max_value=len(files))
I=0
line=4861.325 
step=0.015
ll,rl=4850,4870
c1,c2=60,-60
ll,rl=line-step*int((line-ll)/step),line+step*int((rl-line)/step)
fig,axs=plt.subplots(2)
snr_min=0
for s in files:
    #rv0=RV0(line,s,verbose=False,line=line)
    try:
        SM=sp.spec(s,line=line)
        print(SM)
        rv0=SM.fitline(line)['RV_kms']
        SM.resamp(step,lwl=ll,rwl=rl)
        wv,fl=SM.wave,SM.flux
        pol=np.polyfit(np.concatenate((wv[:c1],wv[c2:])),np.concatenate((fl[:c1],fl[c2:])),1)
        corr=fl/np.polyval(pol,wv)   
        if 1/np.nanstd(np.concatenate((corr[:c1],corr[c2:])))>snr_min: # We avoid data with really bad SNR
            #CRR=cosmic_telluric(corr)
            X_int.append(wv),Y_int.append(corr),HJD.append(SM.hjd),RV.append(rv0)
            axs[0].plot(wv,corr)
        SM.cosmic
        wv,fl=SM.wave,SM.flux
        pol=np.polyfit(np.concatenate((wv[:c1],wv[c2:])),np.concatenate((fl[:c1],fl[c2:])),1)
        corr=fl/np.polyval(pol,wv)   
        if 1/np.nanstd(np.concatenate((corr[:c1],corr[c2:])))>snr_min: # We avoid data with really bad SNR
            #CRR=cosmic_telluric(corr)
            X_int.append(wv),Y_int.append(corr),HJD.append(SM.hjd),RV.append(rv0)
            axs[1].plot(wv,corr)            
    except: 
        pass
    I=I+1
    bar.update(I)






































