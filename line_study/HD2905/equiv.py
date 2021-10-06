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

def fm_fm(line2,line1='4552',star='HD2905'):
    correspondence_FIES,correspondence_MERCATOR,correspondence_SONG=[],[],[]
    lin1=reading_summary(star,line1)
    lin2=reading_summary(star,line2)
    for fil in lin1.Original_file:
        try:
            i=np.where(np.array(lin2.Original_file)==fil)[0][0]
            o=np.where(np.array(lin1.Original_file)==fil)[0][0]
            if  lin1.Instrument[o]=="'FIES'":
                correspondence_FIES.append([lin1.first_moment[o],lin2.first_moment[i]])
            elif lin1.Instrument[o]=="'MERCATOR'":
                correspondence_MERCATOR.append([lin1.first_moment[o],lin2.first_moment[i]])
            elif lin1.Instrument[o]=="'SONG'":
                correspondence_SONG.append([lin1.first_moment[o],lin2.first_moment[i]])
        except:
            pass
    try:
        correspondence_SONG=np.vstack(correspondence_SONG)
    except:
        pass
    try:
        correspondence_MERCATOR=np.vstack(correspondence_MERCATOR)
    except:
        pass        
    try:
        correspondence_FIES=np.vstack(correspondence_FIES)
    except:
        pass
    return correspondence_MERCATOR,correspondence_FIES,correspondence_SONG

# In[prueba_saltable ]

lin=os.listdir('../stars/HD2905/4552')
lines=[]
for l in lin: 
    if l.startswith('HD2905'):
        lines.append(l)
WV_M,V_M,F_M=[],[],[]
WV_F,V_F,F_F=[],[],[]
MERCATOR,FIES=[],[]
for l in lines:
    hdu=open('../stars/HD2905/4552/'+l,'r')
    wave,v,flux=[],[],[]
    for j in hdu:
        dat=j.split()
        wave.append(float(dat[0])),v.append(float(dat[1])),flux.append(float(dat[2]))
    if 'MERCATOR' in l:
        WV_M.append(wave),V_M.append(v),F_M.append(flux)
        l=l.replace('_MERCATOR_4552p622.txt','')
        l=l.replace('HD2905','')
        MERCATOR.append(float(l.replace('p','.')))
    if 'FIES' in l:
        WV_F.append(wave),V_F.append(v),F_F.append(flux)
        l=l.replace('_FIES_4552p622.txt','')
        l=l.replace('HD2905','')
        FIES.append(float(l.replace('p','.')))
    hdu.close()

plt.close('all')
fig,ax=plt.subplots()
pt=[ax.plot(WV_M[i],F_M[i],'r') for i in range(len(WV_M))]
pt=[ax.plot(WV_F[i],F_F[i],'g') for i in range(len(WV_F))]
ax.grid()

WV_M,V_M,F_M=np.vstack(WV_M),np.vstack(V_M),np.vstack(F_M)
WV_F,V_F,F_F=np.vstack(WV_F),np.vstack(V_F),np.vstack(F_F)

MERCATOR_mean=np.mean(F_M,0)
fig1,ax1=plt.subplots()
pt=[ax1.plot(WV_M[i],(F_M[i]/MERCATOR_mean)) for i in range(len(WV_M))]
ax1.grid()

#plt.imshow(MERCATOR,F_M/MERCATOR_mean)

#Z=sorted(zip(MERCATOR,m1_def,mini))
#BJD_list,m1_def,mini=[Z[i][0] for i in range(0,len(BJD_list))],[Z[i][1] for i in range(0,len(BJD_list))],[Z[i][2] for i in range(0,len(BJD_list))]

data=reading_summary('HD2905',folder='4552')
fig2,ax2=plt.subplots()
pt=ax2.plot(np.array(data.HJD)-2450000,data.first_moment,'.')
ax2.grid()

# In[mean_computation]
runcell(0, '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/equiv.py')
def mean_line(star,folder,line,instrument):
    '''
    Parameters
    ----------
    star : str,
        DESCRIPTION.
    folder : str,
        DESCRIPTION.
    line : float,
        DESCRIPTION.
    instrument : str, list or tuple,
        DESCRIPTION.

    Returns
    -------
    WV : array,
        DESCRIPTION.
    V : array,
        DESCRIPTION.
    F : array,
        DESCRIPTION.
    HJD : array,
        DESCRIPTION.
    flux_mean : array,
        DESCRIPTION.

    '''
    path='../stars/'+star+'/'+str(folder)
    lin=os.listdir(path)
    lines=[]
    for l in lin: 
        if l.startswith(star):
            lines.append(l)
        WV,V,F,HJD=[],[],[],[]
    for l in lines:
        hdu=open(path+'/'+l,'r')
        wave,v,flux=[],[],[]
        for j in hdu:
            dat=j.split()
            wave.append(float(dat[0])),v.append(float(dat[1])),flux.append(float(dat[2]))
        if type(instrument)==str:
            if instrument in l:
                str_line=str(line).replace('.','p')
                l=l.replace('_'+instrument+'_'+str_line+'.txt','')
                l=l.replace(star,'')
                HJD.append(float(l.replace('p','.')))
                WV.append(wave),V.append(v),F.append(flux)
        else:
            for ins in instrument:
                if ins in l:
                    str_line=str(line).replace('.','p')
                    l=l.replace('_'+ins+'_'+str_line+'.txt','')
                    l=l.replace(star,'')
                    HJD.append(float(l.replace('p','.')))
                    WV.append(wave),V.append(v),F.append(flux)
        hdu.close()
    WV,V,F,HJD=np.vstack(WV),np.vstack(V),np.vstack(F),np.vstack(HJD)
    flux_mean=np.mean(F,0)
    return WV,V,F,HJD,flux_mean

flux_M,date_M,wave_M,v_M=[],[],[],[]
flux_F,date_F,wave_F,v_F=[],[],[],[]
flux_S,date_S,wave_S,v_S=[],[],[],[]
H_beta=4861.32
WV_M,V_M,F_M,HJD_M,MERCATOR_mean=mean_line('HD2905','H_beta_4861',H_beta,'MERCATOR')
WV_F,V_F,F_F,HJD_F,FIES_mean=mean_line('HD2905','H_beta_4861',H_beta,'FIES')
WV_S,V_S,F_S,HJD_S,SONG_mean=mean_line('HD2905','H_beta_4861',H_beta,'SONG')
fig1,ax1=plt.subplots()
pt=ax1.plot(HJD_M,F_M,'r.')
pt=ax1.plot(HJD_F,F_F,'b.')
pt=ax1.plot(HJD_S,F_S,'g.')


for i in range(len(HJD_M)):
    if HJD_M[i]>2457942 and HJD_M[i]<2458215:
        flux_M.append(F_M[i]/MERCATOR_mean),date_M.append(HJD_M[i])        
        wave_M.append(WV_M[i]),v_M.append(V_M[i])

for i in range(len(HJD_F)):
    if HJD_F[i]>2457942 and HJD_F[i]<2458215:
        flux_F.append(F_F[i]/FIES_mean),date_F.append(HJD_F[i])        
        wave_F.append(WV_F[i]),v_F.append(V_F[i])
        
for i in range(len(HJD_S)):
   if HJD_S[i]>2457942 and HJD_S[i]<2458215:
       flux_S.append(F_S[i]/SONG_mean),date_S.append(HJD_S[i])        
       wave_S.append(WV_S[i]),v_S.append(V_S[i])      
#pt=[plt.plot(WV_M[0],flux[i]) for i in range(len(flux))]
if np.shape(flux_M)!=(0,):
    fecha_M=np.zeros(np.shape(flux_M))
    for i in range(len(fecha_M[0,:])):
        for j in range(len(fecha_M[:,0])): 
            fecha_M[j,i]=date_M[j]
    fecha_S=np.zeros(np.shape(flux_S))
if np.shape(flux_S)!=(0,):
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
pt = plt.scatter(wave_S,fecha_S,c=flux_S,s=0.5,cmap=cm,vmax=1.05)
pt = plt.scatter(wave_M,fecha_M,c=flux_M,s=0.5,cmap=cm,vmax=1.05)
plt.show()
plt.colorbar()
plt.title(r'Residues, $\lambda_0=$'+str(4552.662))
plt.ylabel('HJD [days]')
plt.xlabel(r'$\lambda(\AA)$')
if len(wave_M[0])<len(fecha_M[0]):
    plt.xlim(wave_M[0][0],wave_M[0][-1])
else:
    plt.xlim(wave_S[0][0],wave_S[0][-1])


# In[GIF_H_beta]
gif=False
if gif==True:
    Z=sorted(zip(date_S,wave_S,F_S,v_S))
    date,wave,flux=[Z[i][0] for i in range(0,len(date_S))],[Z[i][1] for i in range(0,len(date_S))],[Z[i][2] for i in range(0,len(date_S))]
    velocity=[Z[i][3] for i in range(0,len(date_S))]
    SM=reading_summary('HD2905','H_beta_4861')
    FM,JD=[],[]
    for i in range(len(np.array(SM.HJD))):
        if SM.HJD[i]>2457942 and SM.HJD[i]<2458215:
            JD.append(SM.HJD[i]),FM.append(SM.first_moment[i])
    import sys
    from matplotlib.animation import FuncAnimation
    #wave,flux=np.concatenate((wave_M,wave_S)),np.concatenate((flux_M,flux_S))
    #fig, ax = plt.subplots(21)
    try:
        plt.close(fig)
    except:
        pass
    fig,axes= plt.subplots(2,2)#,sharex=True)
    fig.subplots_adjust(right=0.85)
    #fig.set_tight_layout(True)
    ax1, ax3=axes[0]
    ax2, ax4=axes[1]
    ax2.set_xlabel(r'Velocity [km/s]')
    ax1.set_ylabel(r'Normalized flux')
    ax2.set_ylabel('HJD [days]-2450000')
    cm = plt.cm.get_cmap('RdYlBu')
    im = ax2.scatter(v_S,fecha_S-2.45e6,c=flux_S,s=0.5,cmap=cm,vmax=1.05)
    #ax2.axhline(y=sorted(fecha_S[:,0])[0]-2.45e6,c='r')
    #plt.colorbar(im,ax=ax2)
    # Query the figure's on-screen size and DPI. Note that when saving the figure to
    # a file, we need to provide a DPI for that separately.
    print('fig size: {0} DPI, size in inches {1}'.format(
        fig.get_dpi(), fig.get_size_inches()))
    #ax3.axvline(x=sorted(fecha_S[:,0])[0]-2.45e6,c='r')
    # Plot a scatter that persists (isn't redrawn) and the initial line.
    x = velocity[0] 
    #ax1.scatter(x, x + np.random.normal(0, 3.0, len(x)))
    line, = ax1.plot(x, flux[0], 'r-', linewidth=2)
    line2 = ax2.axhline(y=sorted(fecha_S[:,0])[0]-2.45e6,c='r')
    #ax1.sharex(ax2)
    ax1.grid()
    ax1.set_ylim(0.7,1.2)
    ax1.plot(x,SONG_mean,'grey')
    #ax1.plot(x,np.mean((flux),0),'g')
    ax2.set_xlim(min(velocity[0]),max(velocity[0]))
    ax1.set_xlim(min(velocity[0]),max(velocity[0]))
    # Para compartir eje
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    
    ax4.plot(np.array(JD)-2.45e6,FM,'.')
    line3 = ax4.axvline(x=np.array(sorted(JD)[0])-2.45e6,c='r')
    ax3.grid()
    ax4.set_xlabel('HJD [days]-2450000')
    ax3.plot(np.array(SM.HJD)-2.45e6,SM.first_moment,'.')
    ax3.axvline(x=sorted(fecha_S[:,0])[0]-2.45e6,c='g')
    ax3.axvline(x=sorted(fecha_S[:,0])[-1]-2.45e6,c='g')
    ax4.grid()
    ax3.set_ylabel(r'$\langle v\rangle$ [km/s]',labelpad=-5),ax4.set_ylabel(r'$\langle v\rangle$ [km/s]',labelpad=-5)
    #cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.95)
    fig.tight_layout(w_pad=6)
    cb_ax = fig.add_axes([0.45, 0.097, 0.02, 0.395]) #[left, bottom, width, height] 
    cbar = fig.colorbar(im, cax=cb_ax)
    #cbar.set_ticks(np.arange(0, 2, 0.5))
    #cbar.set_ticklabels(['low', 'medium', 'high'])
    
    plt.show()
    #ax1.set_xscale('log')  
    #ax1.set_yscale('log')
    #ax1.set_ylim(1e-15,1e-1)
    bar = progressbar.ProgressBar(max_value=np.shape(flux)[0])
    def update(i):
        #label = 'shell {0}'.format(i+Kgen_8000[0][0])
        title = r'H$\beta$, '+str(date[i][0].round(3))
        #print(label)
        # Update the line and the ax1es (with a new xlabel). Return a tuple of
        # "artists" that have to be redrawn for this frame.
        line.set_ydata(flux[i])
        line2.set_ydata(sorted(fecha_S[:,0])[i]-2.45e6)
        line3.set_xdata(np.array(sorted(JD)[i])-2.45e6)
        fig.suptitle(title)
        fig.subplots_adjust(top=0.95)
        bar.update(i)
        return line, ax1
    if __name__ == '__main__':
        # FuncAnimation will call the 'update' function for each frame; here
        # animating over 10 frames, with an interval of 200ms between frames.
        anim = FuncAnimation(fig, update, frames=np.shape(flux)[0], interval=200)
        if len(sys.argv) > 1 and sys.argv[1] == 'save':
            anim.save('line.gif', dpi=80, writer='imagemagick')
        else:
            # plt.show() will just loop the animation forever.
            plt.show()    
        anim.save('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/GIFS/H_beta_HD2905.gif', writer='imagemagick', fps=3)

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

axs[0,1].plot(HJD_SiIII-2.45e6,ZM_SiIII-np.mean(ZM_SiIII),'s',label='Si III')
axs[0,1].plot(HJD_OII-2.45e6,ZM_OII-np.mean(ZM_OII),'^',label='O II')
axs[0,1].plot(HJD_Halpha-2.45e6,ZM_Halpha-np.mean(ZM_Halpha),'.',label=r'H$\alpha$')
axs[0,1].plot(HJD_Hbeta-2.45e6,ZM_Hbeta-np.mean(ZM_Hbeta),'.',label=r'H$\beta$')
axs[0,1].plot(HJD_HeI-2.45e6,ZM_HeI-np.mean(ZM_HeI),'.',label='He I')
axs[0,1].grid(),axs[0,1].legend()#,axs[0,1].set_ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
axs[0,1].set_xlim(6640,6650)#,axs[0,1].set_ylim(-75,25)


axs[0,2].plot(HJD_SiIII-2.45e6,ZM_SiIII-np.mean(ZM_SiIII),'s',label='Si III')
axs[0,2].plot(HJD_OII-2.45e6,ZM_OII-np.mean(ZM_OII),'^',label='O II')
axs[0,2].plot(HJD_Halpha-2.45e6,ZM_Halpha-np.mean(ZM_Halpha),'.',label=r'H$\alpha$')
axs[0,2].plot(HJD_Hbeta-2.45e6,ZM_Hbeta-np.mean(ZM_Hbeta),'.',label=r'H$\beta$')
axs[0,2].plot(HJD_HeI-2.45e6,ZM_HeI-np.mean(ZM_HeI),'.',label='He I')
axs[0,2].grid(),axs[0,2].legend()#,axs[0,2].set_xlabel('HJD [days]-2450000')#,axs[0,2].set_ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
axs[0,2].set_xlim(7575,7620)#,axs[0,2].set_ylim(-50,50)


SiIII=reading_summary('../stars/HD2905/4552/00_HD2905_4552.txt')
SiIV=reading_summary('../stars/HD2905/4116/00_HD2905_4116.txt')
OII=reading_summary('../stars/HD2905/4661/00_HD2905_4661.txt')

axs[1,0].plot(np.vstack(SiIII.HJD)-2.45e6,SiIII.first_moment,'s',label='Si III')
axs[1,0].plot(np.vstack(SiIV.HJD)-2.45e6,SiIV.first_moment,'.',label='Si IV')
axs[1,0].plot(np.vstack(OII.HJD)-2.45e6,OII.first_moment,'.',label='O II')
axs[1,0].legend(),axs[1,0].grid(),axs[1,0].set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$ [km/s]')

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



# In[Telescope_effect_SiIII]


SM=reading_summary('../stars/HD2905/4574/00_HD2905_4574.txt')
m_SiIII_fm,f_SiIII_fm,s_SiIII_fm=[],[],[]
m_HJD,f_HJD,s_HJD=[],[],[]
for i in range(len(SM.first_moment)):
    if SM.Instrument[i]=="'MERCATOR'":
        m_SiIII_fm.append(SM.first_moment[i]),m_HJD.append(SM.HJD[i])
    elif SM.Instrument[i]=="'FIES'":
        f_SiIII_fm.append(SM.first_moment[i]),f_HJD.append(SM.HJD[i])
    elif SM.Instrument[i]=="'SONG'":
        s_SiIII_fm.append(SM.first_moment[i]),s_HJD.append(SM.HJD[i])
    else:
        print('Unknown instrument')

plt.plot(np.vstack(m_HJD)-2.45e6,m_SiIII_fm,'xr',label='HERMES')        
plt.plot(np.vstack(f_HJD)-2.45e6,f_SiIII_fm,'xb',label='FIES')            
plt.plot(np.vstack(s_HJD)-2.45e6,s_SiIII_fm,'xg',label='SONG')

plt.plot(np.vstack(m_HJD)-2.45e6,np.vstack(m_SiIII_fm)-np.mean(m_SiIII_fm),'.r',label='Rescaled HERMES')        
plt.plot(np.vstack(f_HJD)-2.45e6,np.vstack(f_SiIII_fm)-np.mean(f_SiIII_fm),'.b',label='Rescaled FIES')            
plt.plot(np.vstack(s_HJD)-2.45e6,np.vstack(s_SiIII_fm)-np.mean(s_SiIII_fm),'.g',label='Rescaled SONG')
plt.grid(),plt.legend()

plt.ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
plt.xlabel('HJD [days]-2450000')














