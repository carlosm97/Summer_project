#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 14:26:58 2021

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
# Let's study 2 or 3 lines each star. 
def general(star,folder_fm='4567',ZM=True):
    folders=['H_alpha_6562','H_alpha']
    lines=[6562.80,6562.80]
    if ZM==True:
        for i in range(len(lines)):
            try:
                if os.path.exists('../stars/'+star+'/'+folders[i]+'/00_zero_moment.txt')==True: 
                    os.remove('../stars/'+star+'/'+folders[i]+'/00_zero_moment.txt')
                SM=reading_summary('../stars/'+star+'/'+folders[i]+'/00_'+star+'_'+folders[i]+'.txt',verbose=False)
                files=os.listdir('../stars/'+star+'/'+folders[i])
                f=open('../stars/'+star+'/'+folders[i]+'/00_zero_moment.txt', 'a')
                for fil in range(len(SM.line_file)):
                    wv,v,fl=reading_line_file('../stars/'+star+'/'+folders[i]+'/'+SM.line_file[fil][1:-1])
                    l1,l2=fmr.find_nearest(wv,SM.lam_min[fil])[1],fmr.find_nearest(wv,SM.lam_max[fil])[1]
                    info=SM.Original_file[fil][2:-1],SM.line_file[fil][1:-1],SM.HJD[fil],zero_moment(wv[l1:l2],fl[l1:l2],lines[i])
                    f.writelines(str(info))
                    f.writelines('\n')
                f.close()
            except: print('Passing ',folders[i])
    try:    
        OF,LF,HJD,ZM=zero_reading(star,'H_alpha_6562')
    except:
        OF,LF,HJD,ZM=zero_reading(star,'H_alpha')
    FM=reading_summary(star,folder_fm)
    
    HJD_ph,mag=photometry(star)
    fig,axs=plt.subplots(3,sharex=True)
    axs[0].plot(HJD-2.45e6,ZM,'.')
    axs[1].plot(np.array(FM.HJD)-2.45e6,FM.first_moment,'.')
    axs[2].plot(HJD_ph,mag,'.',markersize=0.5)
    axs[0].grid(),axs[1].grid(),axs[2].grid()
    axs[2].set_xlabel(r'HJD[days]-2450000')
    axs[0].set_ylabel(r'$\langle v^0\rangle, H\alpha$')
    axs[1].set_ylabel(r'$\langle v\rangle, SiIII(4567.84)$')
    axs[2].set_ylabel('mag[mmag]')
'''
stars=[]
for s in os.listdir('../stars'):
    if s.startswith('HD'):
        stars.append(s)
'''


7
plt.rcParams["figure.figsize"] = (15,12)
def line_summary(star):
    plt.rcParams["figure.figsize"] = (15,12)
    #star='HD2905'
    spectrum=findstar(star)
    files=[ins for ins in spectrum if 'FH' in ins or 'SONG' in ins]
    plt.close('all')
    c= 299792.458
    fig, axs=plt.subplots(3,3)
    step=0.015
    
    def line_plot(s,ln,l1,l2):
        try: SM=sp.spec(s,line=ln)
        except: return None
        SM.resamp(step,lwl=l1,rwl=l2)
        wv,fl=SM.wave,SM.flux
        v=(np.array(wv)-ln)*c/ln
        if 'SONG' in s:
            pol=np.polyfit(np.concatenate((wv[:50],wv[50:])),np.concatenate((fl[:50],fl[50:])),1)
            fl=fl/np.polyval(pol,wv) 
            if 1/np.nanstd(np.concatenate((fl[:50],fl[50:])))<15:
                return None
        return np.array([v,fl]) 
    H_alpha_l,H_beta_l,HeI_5875_l,SiIV_l,SiIII_l,SiII_l,OII_l,OIII_l,HeI_5015_l=6562.80,4861.325,5875.62,\
    4116.103,4567.840,6347.11,4661.6324,5592.252,5015.678
    bar = progressbar.ProgressBar(max_value=len(files))
    I=0
    H_alpha_v,H_beta_v,HeI_5875_v,SiIV_v,SiIII_v,SiII_v,OII_v,OIII_v,HeI_5015_v=[],[],[],[],[],[],[],[],[]  
    H_alpha_f,H_beta_f,HeI_5875_f,SiIV_f,SiIII_f,SiII_f,OII_f,OIII_f,HeI_5015_f=[],[],[],[],[],[],[],[],[]    
    ln_Ha,ln_Hb,ln_5875,ln_5015,ln_SiII,ln_SiIII,ln_SiIV,ln_OII,ln_OIII=0,0,0,0,0,0,0,0,0
    for s in files:
        lin=H_alpha_l
        line=line_plot(s,lin,lin-10,lin+10)
        if np.any(line)!=None:
            H_alpha_v.append(line[0]),H_alpha_f.append(line[1])
            ln_Ha+=1

        lin=H_beta_l
        line=line_plot(s,lin,lin-10,lin+10)   
        if np.any(line)!=None: 
            pol=np.polyfit(np.concatenate((line[0][:60],line[0][-60:])),np.concatenate((line[1][:60],line[1][-60:])),1)
            line[1]=line[1]/np.polyval(pol,line[0])  
            H_beta_v.append(line[0]),H_beta_f.append(line[1])
            ln_Hb+=1

        lin=HeI_5875_l
        line=line_plot(s,lin,lin-10,lin+10) 
        if np.any(line)!=None:
            HeI_5875_v.append(line[0]),HeI_5875_f.append(line[1])
            ln_5875+=1
        
        lin=SiIV_l
        line=line_plot(s,lin,lin-5,lin+5) 
        if np.any(line)!=None:
            SiIV_v.append(line[0]),SiIV_f.append(line[1])
            ln_SiIV+=1
        
        lin=SiIII_l
        line=line_plot(s,lin,lin-5,lin+5)
        if np.any(line)!=None:
            pol=np.polyfit(np.concatenate((line[0][:60],line[0][-60:])),np.concatenate((line[1][:60],line[1][-60:])),1)
            line[1]=line[1]/np.polyval(pol,line[0])  
            SiIII_v.append(line[0]),SiIII_f.append(line[1])
            ln_SiIII+=1
        
        if 'SONG' not in s:
            lin=SiII_l
            line=line_plot(s,lin,lin-5,lin+5)  
            if np.any(line)!=None:
                SiII_v.append(line[0]),SiII_f.append(line[1])
                ln_SiII+=1
            
            lin=OIII_l
            line=line_plot(s,lin,lin-5,lin+5) 
            if np.any(line)!=None:
                OIII_v.append(line[0]),OIII_f.append(line[1])
                ln_OIII+=1
        
        lin=OII_l
        line=line_plot(s,lin,lin-3.5,lin+3.5)  
        if np.any(line)!=None:
            pol=np.polyfit(np.concatenate((line[0][:60],line[0][-60:])),np.concatenate((line[1][:60],line[1][-60:])),1)
            line[1]=line[1]/np.polyval(pol,line[0])  
            OII_v.append(line[0]),OII_f.append(line[1])
            ln_OII+=1

        lin=HeI_5015_l
        line=line_plot(s,lin,lin-4,lin+4) 
        if np.any(line)!=None:
            pol=np.polyfit(np.concatenate((line[0][:60],line[0][-60:])),np.concatenate((line[1][:60],line[1][-60:])),1)
            line[1]=line[1]/np.polyval(pol,line[0])  
            HeI_5015_v.append(line[0]),HeI_5015_f.append(line[1])
            ln_5015+=1
        
        I=I+1
        bar.update(I)
    H_alpha_v_m,H_alpha_f_m=np.mean(H_alpha_v,0),np.mean(H_alpha_f,0)
    H_beta_v_m,H_beta_f_m=np.mean(H_beta_v,0),np.mean(H_beta_f,0)
    HeI_5875_v_m,HeI_5875_f_m=np.mean(HeI_5875_v,0),np.mean(HeI_5875_f,0)
    SiIV_v_m,SiIV_f_m=np.mean(SiIV_v,0),np.mean(SiIV_f,0)
    SiIII_v_m,SiIII_f_m=np.mean(SiIII_v,0),np.mean(SiIII_f,0)
    SiII_v_m,SiII_f_m=np.mean(SiII_v,0),np.mean(SiII_f,0)
    OII_v_m,OII_f_m=np.mean(OII_v,0),np.mean(OII_f,0)
    OIII_v_m,OIII_f_m=np.mean(OIII_v,0),np.mean(OIII_f,0)
    HeI_5015_v_m,HeI_5015_f_m=np.mean(HeI_5015_v,0),np.mean(HeI_5015_f,0)
    if ln_Ha>=50: rnd=np.random.randint(0,ln_Ha,50)
    else: rnd=range(ln_Ha)
    [axs[0,0].plot(H_alpha_v[rn],H_alpha_f[rn],alpha=0.5) for rn in rnd]

    if ln_Hb>=50: rnd=np.random.randint(0,ln_Hb,50)
    else: rnd=range(ln_Hb)
    [axs[0,1].plot(H_beta_v[rn],H_beta_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_5875>=50: rnd=np.random.randint(0,ln_5875,50)
    else: rnd=range(ln_5875)
    [axs[0,2].plot(HeI_5875_v[rn],HeI_5875_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_SiIV>=50: rnd=np.random.randint(0,ln_SiIV,50)
    else: rnd=range(ln_SiIV)
    [axs[1,0].plot(SiIV_v[rn],SiIV_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_SiIII>=50: rnd=np.random.randint(0,ln_SiIII,50)
    else: rnd=range(ln_SiIII)
    [axs[1,1].plot(SiIII_v[rn],SiIII_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_SiII>=50: rnd=np.random.randint(0,ln_SiII,50)
    else: rnd=range(ln_SiII)
    [axs[1,2].plot(SiII_v[rn],SiII_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_OII>=50: rnd=np.random.randint(0,ln_OII,50)
    else: rnd=range(ln_OII)
    [axs[2,0].plot(OII_v[rn],OII_f[rn],alpha=0.5) for rn in rnd]
    
    if ln_OIII>=50: rnd=np.random.randint(0,ln_OIII,50)
    else: rnd=range(ln_OIII)
    [axs[2,1].plot(OIII_v[rn],OIII_f[rn],alpha=0.5) for rn in rnd]

    if ln_5015>=50: rnd=np.random.randint(0,ln_5015,50)
    else: rnd=range(ln_5015)
    [axs[2,2].plot(HeI_5015_v[rn],HeI_5015_f[rn]) for rn in rnd]
    axs[0,0].set_ylim(min(H_alpha_f_m)-0.5*(max(H_alpha_f_m)-min(H_alpha_f_m)),max(H_alpha_f_m)+0.5*(max(H_alpha_f_m)-min(H_alpha_f_m)))
    axs[0,1].set_ylim(min(H_beta_f_m)-0.5*(max(H_beta_f_m)-min(H_beta_f_m)),max(H_beta_f_m)+0.5*(max(H_beta_f_m)-min(H_beta_f_m)))
    axs[0,2].set_ylim(min(HeI_5875_f_m)-0.5*(max(HeI_5875_f_m)-min(HeI_5875_f_m)),max(HeI_5875_f_m)+0.5*(max(HeI_5875_f_m)-min(HeI_5875_f_m)))
    axs[1,0].set_ylim(min(SiIV_f_m)-0.5*(max(SiIV_f_m)-min(SiIV_f_m)),max(SiIV_f_m)+0.5*(max(SiIV_f_m)-min(SiIV_f_m)))
    axs[1,1].set_ylim(min(SiIII_f_m)-0.5*(max(SiIII_f_m)-min(SiIII_f_m)),max(SiIII_f_m)+0.5*(max(SiIII_f_m)-min(SiIII_f_m)))
    axs[1,2].set_ylim(min(SiII_f_m)-0.5*(max(SiII_f_m)-min(SiII_f_m)),max(SiII_f_m)+0.5*(max(SiII_f_m)-min(SiII_f_m)))
    axs[2,0].set_ylim(min(OII_f_m)-0.5*(max(OII_f_m)-min(OII_f_m)),max(OII_f_m)+0.5*(max(OII_f_m)-min(OII_f_m)))
    axs[2,1].set_ylim(min(OIII_f_m)-0.5*(max(OIII_f_m)-min(OIII_f_m)),max(OIII_f_m)+0.5*(max(OIII_f_m)-min(OIII_f_m)))
    axs[2,2].set_ylim(min(HeI_5015_f_m)-0.5*(max(HeI_5015_f_m)-min(HeI_5015_f_m)),max(HeI_5015_f_m)+0.5*(max(HeI_5015_f_m)-min(HeI_5015_f_m)))
   
    axs[0,0].plot(H_alpha_v_m,H_alpha_f_m,'k',markersize=0.5)
    axs[0,1].plot(H_beta_v_m,H_beta_f_m,'k',markersize=0.5)
    axs[0,2].plot(HeI_5875_v_m,HeI_5875_f_m,'k',markersize=0.5)
    axs[1,0].plot(SiIV_v_m,SiIV_f_m,'k',markersize=0.5)
    axs[1,1].plot(SiIII_v_m,SiIII_f_m,'k',markersize=0.5)
    axs[1,2].plot(SiII_v_m,SiII_f_m,'k',markersize=0.5)
    axs[2,0].plot(OII_v_m,OII_f_m,'k',markersize=0.5)
    axs[2,1].plot(OIII_v_m,OIII_f_m,'k',markersize=0.5)
    axs[2,2].plot(HeI_5015_v_m,HeI_5015_f_m,'k',markersize=0.5)
   
    axs[0,0].grid(),axs[0,1].grid(),axs[0,2].grid(),axs[1,0].grid(),axs[1,1].grid(),
    axs[1,2].grid(),axs[2,0].grid(),axs[2,1].grid(),axs[2,2].grid()
    axs[0,0].set_title(r'$H\alpha$, '+str(round(H_alpha_l,3)))
    axs[0,1].set_title(r'$H\beta$, '+str(round(H_beta_l,3)))
    axs[0,2].set_title(r'$HeI$, '+str(round(HeI_5875_l,3)))
    axs[1,0].set_title(r'$SiIV$, '+str(round(SiIV_l,3)))
    axs[1,1].set_title(r'$SiIII$, '+str(round(SiIII_l,3)))
    axs[1,2].set_title(r'$SiII$, '+str(round(SiII_l,3)))
    axs[2,0].set_title(r'$OII$, '+str(round(OII_l,3)))
    axs[2,1].set_title(r'$OIII$, '+str(round(OIII_l,3)))
    axs[2,2].set_title(r'$HeI$, '+str(round(HeI_5015_l,3)))      
    axs[1,0].set_ylabel('Normalized flux'),axs[2,1].set_xlabel('Velocity [km/s]')
    plt.suptitle(star,fontsize=15)
    #axs[0,0].set_title(r'$H\alpha$'),axs[0,1].set_title(r'$H\beta$'),axs[0,2].set_title(r'$HeI$ 5875')
    #axs[1,0].set_title(r'$SiIV$'),axs[1,1].set_title(r'$SiIII$'),axs[1,2].set_title(r'$SiII$')
    #axs[2,0].set_title(r'$OII$'),axs[2,1].set_title(r'$OIII$'),axs[2,2].set_title(r'$HeI$ 5015')
  
    path='../stars/SUMMARY'
    if os.path.exists(path)==False:
        os.mkdir(path)
    plt.savefig(path+'/'+star+'_lines_summary.png')
    return None
# In[HD37128]
star='HD37128' #B0Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.85,\
              lim2=1,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'SONG',H_alpha,H_alpha-10,H_alpha+10,0.015,snr_min=15,\
              plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])

SiIII=4567.840
SiIV=4116.103
# We have both.
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,telluric=True,sel='vel',sig=5)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,telluric=True,sel='vel',sig=5)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,plt_mean=True,telluric=True,\
              automatic_selection=False,ll1=-120,ll2=200,sel='vel',sig=5)

lin=SiIV
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=0.9,lim2=1.1,normalization=False,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,automatic_selection=False,\
              ll1=-140,ll2=140,sel='vel',normalization=False,plt_mean=True,telluric=True,sig=3)

general(star)

# In[HD213087]
star='HD213087' #B0.5Ib
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
       
H_alpha=6562.80
SiIII=4567.840
SiIV=4116.103

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.5,\
              lim2=0.55,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])


# We have both.
lin=SiIII
st=line_study(star,'FIES',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=False,\
              sel='vel',ll1=-160,ll2=160,sig=5)
st=line_study(star,'MERCATOR',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=False,\
              sel='vel',ll1=-160,ll2=160,sig=5)

lin=SiIV
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,normalization=False,plt_mean=True,\
              automatic_selection=False,sel='vel',ll1=-150,ll2=105,sig=3)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,automatic_selection=False,\
              ll1=-150,ll2=105,sel='vel',normalization=False,plt_mean=True,sig=3)

general(star)
# In[HD36486]
star='HD36486' #B0III+O9V
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,snr_min=15,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4])

OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True)
st=line_study(star,'SONG',lin,5590,lin+5,0.015,snr_min=15,plt_mean=True,telluric=True)

general(star,'5592')

# In[HD21291]
star='HD21291' #B9Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.5,lim2=0.6,\
              plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,snr_min=15,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)

SiII=6347.11
lin=SiII
prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True,telluric=True)
#st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,snr_min=15,plt_mean=True,telluric=True) There is no line in SONG

SiIII=4552.622
lin=SiIII
st=line_study(star,'FIES',lin,4551,4554,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-60,ll2=50)
st=line_study(star,'MERCATOR',lin,4551,4554,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-60,ll2=50,telluric=True)
    
general(star,'6347')

# In[HD31327]
star='HD31327' #B2.5Ib, SIMULTANEOUS
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.4,lim2=0.6,\
              plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.5,lim2=0.6,\
              plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.8,snr_min=15,\
              plt_mean=True, telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)

    
SiII=6347.11
SiIII=4567.840
SiIV=4116.103

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.5,\
              lim2=0.55,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])


# We have both.
lin=SiIII
st=line_study(star,'FIES',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=False,\
              sel='vel',ll1=-95,ll2=95,sig=5)
st=line_study(star,'MERCATOR',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=False,\
              sel='vel',ll1=-95,ll2=95,sig=5)
st=line_study(star,'SONG',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=False,\
              sel='vel',ll1=-95,ll2=95,sig=4,telluric=True,snr_min=15)

lin=SiII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,normalization=False,plt_mean=True,\
              automatic_selection=False,sel='vel',ll1=-95,ll2=95,sig=3)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,automatic_selection=False,\
              ll1=-95,ll2=95,sel='vel',normalization=False,plt_mean=True,sig=3)
    
general(star)

# In[HD8065]

star='HD8065' #B9Iab Bibliographic RV=-75.3, SIMULTANEOUS
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

#st=line_study(star,'FIES',H_alpha,H_alpha-13,H_alpha+9.5,0.015,lim1=0.5,lim2=0.6,\
#              plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],save=False,\
#                  normalization=False)  #Useless to do statistics
st=line_study(star,'SONG',H_alpha,H_alpha-11.5,H_alpha+9,0.015,snr_min=15,plt_mean=True,\
              automatic_selection=False,ll1=-400,ll2=200,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'FIES',H_alpha,H_alpha-11.5,H_alpha+9,0.015,snr_min=15,plt_mean=True,\
              automatic_selection=False,ll1=-400,ll2=200,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
FeII=6456.383 
lin=FeII
    
st=line_study(star,'SONG',lin,lin-4,lin,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-125,ll2=-30,sig=2,telluric=True,snr_min=15)    
st=line_study(star,'FIES',lin,lin-4,lin,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-125,ll2=-30,sig=2)        
general(star,'6456')

# In[HD35921]

star='HD35921' # O9II
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,lim1=0.5,\
              lim2=0.55,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4])

    
SiIII=4567.84
lin=SiIII
st=line_study(star,'FIES',lin,lin-4,lin+4,0.015,plt_mean=True,automatic_selection=True,\
              sel='vel',ll1=-95,ll2=95,sig=5)
st=line_study(star,'MERCATOR',lin,lin-4,lin+4,0.015,lim1=0.3,lim2=0.4,plt_mean=True,automatic_selection=True,\
              sel='vel',ll1=-95,ll2=95,sig=5)
general(star)

# In[HD34656]

star='HD34656' # O7.5II
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))    
    
H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-11,H_alpha+9.,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-11,H_alpha+9.,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)

OIII=5592.252 
lin=OIII
st=line_study(star,'FIES',lin,lin-4,lin+4,0.015,plt_mean=True,\
              automatic_selection=False,ll1=-150,ll2=150,sig=5)
st=line_study(star,'MERCATOR',lin,lin-4,lin+4,0.015,plt_mean=True,telluric=True,\
          automatic_selection=False,ll1=-150,ll2=150,sig=5)

general(star,'5592')   
    
# In[HD34085]

star='HD34085' # B8Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5,\
              automatic_selection=False,ll1=-200,ll2=200)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,snr_min=15,\
              lim1=0.35,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)





SiII=6347.11
lin=SiII
#prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True,\
              automatic_selection=False,ll1=-75,ll2=100)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,\
              telluric=True,automatic_selection=False,ll1=-75,ll2=100)
#st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,lim2=1.4,plt_mean=True,\
#              telluric=True,automatic_selection=False,ll1=-50,ll2=75)
    
    
    
SiIII=4567.840
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-40,ll2=80)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-40,ll2=80,telluric=True)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-40,ll2=80,telluric=True,snr_min=15)

general(star,'4567')   

# In[HD198478]

star='HD198478' # B4Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5,\
              lim1=0.35,lim2=0.37,normalization=False)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)
st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.45,0.015,snr_min=15,\
              lim1=0.75,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5,save=False)




    
    
SiIII=4567.840
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.8,lim2=1)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.9,lim2=1)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,plt_mean=True,telluric=True,snr_min=15,sig=5)

general(star,'4567')   


SiII=6347.11
lin=SiII
#prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=0.65,sig=5,plt_mean=True,telluric=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=0.7,sig=5,plt_mean=True,telluric=True)
# In[HD199478]

star='HD199478' # B8Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80
st=line_study(star,'FIES',H_alpha,H_alpha-9.5,H_alpha+9.4,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sel='vel',sig=5,\
              lim1=0.7,lim2=0.8,normalization=False)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-9.5,H_alpha+9.4,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5,lim1=0.6,lim2=0.8,normalization=False)



SiIII=4567.840
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-90,ll2=50)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-90,ll2=50)
general(star,'4567')   


SiII=6347.11
lin=SiII
#prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=0.35,lim2=0.4,sig=5,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=0.6,lim2=0.7,sig=5,plt_mean=True)

# In[HD36512]

star='HD36512' # O9.7V
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80
st=line_study(star,'FIES',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,\
              automatic_selection=False,ll1=-750,ll2=750)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,\
              automatic_selection=False,ll1=-750,ll2=750)
#st=line_study(star,'SONG',H_alpha,H_alpha-9.5,H_alpha+9.5,0.015,snr_min=15,\
#              lim1=0.75,plt_mean=True,telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5)

OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-35,ll2=60)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-35,ll2=60,telluric=True) 
st=line_study(star,'SONG',lin,lin-2,lin+2,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-35,ll2=60,telluric=True,snr_min=15)
general(star,'4567')   


OII=4661.6324
lin=OII

st=line_study(star,'FIES',lin,lin-3,lin+3,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-30,ll2=65)
st=line_study(star,'MERCATOR',lin,lin-3,lin+3,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-30,ll2=65,telluric=True) 
st=line_study(star,'SONG',lin,lin-3,lin+3,0.015,plt_mean=True,automatic_selection=False,\
          ll1=-30,ll2=65,snr_min=15,telluric=True) 
general(star,'5592')   


# In[HD214680]

star='HD214680' # O9V
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80
st=line_study(star,'FIES',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,\
              automatic_selection=False,ll1=-750,ll2=750)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,\
              automatic_selection=False,ll1=-750,ll2=750)

OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-2.5,lin+2.5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-75,ll2=50)
st=line_study(star,'MERCATOR',lin,lin-2.5,lin+2.5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-75,ll2=50,telluric=True) 
st=line_study(star,'SONG',lin,lin-2.5,lin+2.5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-75,ll2=50,telluric=True,snr_min=15)     
general(star,'4567')   

OII=4661.6324
lin=OII

st=line_study(star,'FIES',lin,lin-2,lin+1.4,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-60,ll2=40)
st=line_study(star,'MERCATOR',lin,lin-2,lin+1.4,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-60,ll2=40,telluric=True) 
st=line_study(star,'SONG',lin,lin-2,lin+1.4,0.015,plt_mean=True,automatic_selection=False,\
          ll1=-60,ll2=40,snr_min=15,telluric=True) 
general(star,'5592')   


# In[HD188001]

star='HD188001' # O7.5Iabf
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80


st=line_study(star,'FIES',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-30,H_alpha+30,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],lim1=1.2,lim2=1.4,sig=10)
    
    
    
OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.8,sig=5)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,automatic_selection=False,\
              telluric=True)      
# No OII

general(star,'5592')   



# In[HD47839]

star='HD47839' # O7V+B1.5/2V
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,automatic_selection=False,\
              ll1=-500,ll2=640)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,automatic_selection=False,\
              ll1=-500,ll2=640)



OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.8,sig=5,automatic_selection=False,\
              ll1=-75,ll2=130)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,automatic_selection=False,\
              telluric=True,ll1=-75,ll2=130)      
st=line_study(star,'SONG',lin,5590,lin+5,0.015,plt_mean=True,sig=5,telluric=True,snr_min=15,\
              automatic_selection=False,ll1=-75,ll2=130)       
    
general(star,'5592')     
    
    
# In[HD38771]

star='HD38771' # B0.5Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80

st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,lim1=0.6,lim2=0.8)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,lim1=0.9)
st=line_study(star,'SONG',H_alpha,H_alpha-10,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,\
              snr_min=15,lim1=1.35,lim2=1.6)

SiIII=4567.840
lin=SiIII
#prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,sig=5,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,sig=5,plt_mean=True)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,sig=5,plt_mean=True,snr_min=15,\
              lim1=1.2,lim2=1.3)
general(star,'4567')     

SiIV=4116.103
lin=SiIV
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,sig=5,plt_mean=True,normalization=False,\
              automatic_selection=False,ll1=-150,ll2=150)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,sig=5,plt_mean=True,normalization=False,\
              telluric=True,automatic_selection=False,ll1=-150,ll2=150)
# In[HD190603]

star='HD190603' # B1.5Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 

H_alpha=6562.80    
st=line_study(star,'FIES',H_alpha,H_alpha-10,H_alpha+10,0.015,lim1=1.2,lim2=1.4,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-10,H_alpha+10,0.015,lim1=1.2,lim2=1.4,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10)    
st=line_study(star,'SONG',H_alpha,H_alpha-10,H_alpha+10,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,snr_min=15)    


SiIII=4567.840
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,telluric=True)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,plt_mean=True,snr_min=15,telluric=True)
general(star,'4567')   

SiIV=4116.103
lin=SiIV
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=1.05,sig=5,plt_mean=True,\
              normalization=False,automatic_selection=False,ll1=-50,ll2=110)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=0.6,lim2=0.8,sig=5,plt_mean=True,\
              normalization=False,automatic_selection=False,ll1=-50,ll2=110,telluric=True)
'''
SiII=6347.11
lin=SiII
#prueba_linea(star,'FIES',lin,lin-5,lin+5)
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,lim1=1.05,sig=5,plt_mean=True,telluric=True)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,lim1=0.6,lim2=0.8,sig=5,plt_mean=True,telluric=True)
'''

# In[HD207198]
  
star='HD207198' # O8.5II
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star)))
line_summary(star)


H_alpha=6562.80    
st=line_study(star,'FIES',H_alpha,H_alpha-9.8,H_alpha+10.5,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5,automatic_selection=False,ll1=-400,ll2=400)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-10,H_alpha+10,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=5,automatic_selection=False,ll1=-400,ll2=400)    
st=line_study(star,'SONG',H_alpha,H_alpha-10,H_alpha+10,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],snr_min=15,sig=5,automatic_selection=False,ll1=-400,ll2=400)  

OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-150,ll2=140)
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,automatic_selection=False,\
              ll1=-150,ll2=140,telluric=True) 
general(star,'5592')   


# In[HD206165]

star='HD206165' # B2Ib
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 
line_summary(star)
SB(star)['SP_TYPE'][0]
H_alpha=6562.80


st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False,lim1=0.5,lim2=0.7)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False,lim1=0.5,lim2=0.7)
    
SiIII=4567.84
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,lim1=0.8,sig=5)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,lim1=0.9)     
    
SiII=6347.11
lin=SiII
st=line_study(star,'FIES',lin,lin-3.5,lin+3.5,0.015,plt_mean=True,lim1=0.8,sig=5,\
              automatic_selection=False,ll1=-100,ll2=75)    acc
st=line_study(star,'MERCATOR',lin,lin-3.5,lin+3.5,0.015,plt_mean=True,sig=5,\
              telluric=True,automatic_selection=False,ll1=-100,ll2=75)   
    
general(star,'4567')   

# In[HD30614]

star='HD30614' # O9Ia
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 
line_summary(star)
print(SB(star)['SP_TYPE'][0])
H_alpha=6562.80
st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False,lim1=1.4,lim2=1.6)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False,lim1=1.4,lim2=1.6)
#st=line_study(star,'SONG',H_alpha,H_alpha-10,H_alpha+10,0.015,plt_mean=True,\
#              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,snr_min=15,lim1=1.4,lim2=1.6)    
# Muy ajustado.
OIII=5592.252
lin=OIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,normalization=False)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,normalization=False) 
    
general(star,'5592')
# In[HD209975]
star='HD209975' # O9Ib
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 
line_summary(star)
print(SB(star)['SP_TYPE'][0])
H_alpha=6562.80  
st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False,lim1=1.4,lim2=1.6)
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,target='H_alpha_'+str(H_alpha)[:4],sig=10,normalization=False)
  

SiIII=4567.84
lin=SiIII
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,automatic_selection=False,ll1=-150,ll2=120)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,automatic_selection=False,ll1=-150,ll2=120,telluric=True)



SiIV=4116.103
lin=SiIV   
st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,normalization=False,automatic_selection=False,ll1=-150,ll2=130)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,normalization=False,automatic_selection=False,ll1=-150,ll2=130,telluric=True) 



# In[HD164353]
star='HD164353'	#B5I
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 
line_summary(star)
print(SB(star)['SP_TYPE'][0])
H_alpha=6562.80  
st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,normalization=False,automatic_selection=False,ll1=-600,\
              ll2=600,sig=10,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,normalization=False,automatic_selection=False,ll1=-600,\
              ll2=600,sig=10,target='H_alpha_'+str(H_alpha)[:4])  
    
SiIII=4567.84
lin=SiIII

st=line_study(star,'FIES',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              automatic_selection=False,ll1=-80,ll2=75)    
st=line_study(star,'MERCATOR',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
             snr_min=15,automatic_selection=False,ll1=-80,ll2=75)
st=line_study(star,'SONG',lin,lin-5,lin+5,0.015,plt_mean=True,sig=5,\
              snr_min=15,telluric=True,automatic_selection=False,ll1=-80,ll2=75)    
    
# In[HD24431]   
star='HD24431'
if os.path.exists('../stars/'+star)==False:
    os.mkdir('../stars/'+star)
print(len(findstar(star))) 
line_summary(star)
print(SB(star)['SP_TYPE'][0])
H_alpha=6562.80  
'''
st=line_study(star,'FIES',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
             telluric=True,normalization=False,automatic_selection=False,ll1=-600,\
             ll2=600,sig=10,target='H_alpha_'+str(H_alpha)[:4])
st=line_study(star,'MERCATOR',H_alpha,H_alpha-15,H_alpha+15,0.015,plt_mean=True,\
              telluric=True,normalization=False,automatic_selection=False,ll1=-600,\
              ll2=600,sig=10,target='H_alpha_'+str(H_alpha)[:4]) '''
prueba_linea(star,'FIES',5888,5880,5895)


prueba_linea(star,'FIES',5592,5587,5597)






   