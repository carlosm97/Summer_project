#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 13:48:34 2021

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
from scipy import stats
from matplotlib.patches import Rectangle

def statistical(star,lin_fol):
    SM=reading_summary('../stars/'+star+'/'+lin_fol+'/00_'+star+'_'+lin_fol+'.txt')
    m_SiIII_fm,f_SiIII_fm,s_SiIII_fm=[],[],[]
    m_HJD,f_HJD,s_HJD=[],[],[]
    error_m,error_f,error_s=[],[],[]
    for i in range(len(SM.first_moment)):
        if SM.Instrument[i]=="'MERCATOR'":
            m_SiIII_fm.append(SM.first_moment[i]),m_HJD.append(SM.HJD[i]),error_m.append(SM.error_fm[i])
        elif SM.Instrument[i]=="'FIES'":
            f_SiIII_fm.append(SM.first_moment[i]),f_HJD.append(SM.HJD[i]),error_f.append(SM.error_fm[i])
        elif SM.Instrument[i]=="'SONG'":
            s_SiIII_fm.append(SM.first_moment[i]),s_HJD.append(SM.HJD[i]),error_s.append(SM.error_fm[i])
        else:
            print('Unknown instrument')
    plt.figure()
    #plt.plot(np.vstack(m_HJD)-2.45e6,m_SiIII_fm,'xr',label='HERMES')        
    #plt.plot(np.vstack(f_HJD)-2.45e6,f_SiIII_fm,'xb',label='FIES')            
    #plt.plot(np.vstack(s_HJD)-2.45e6,s_SiIII_fm,'xg',label='SONG')
    
    plt.ylabel(r'$\langle v^0\rangle-\langle v^0_0\rangle$ [km/s]')
    plt.xlabel('HJD [days]-2450000')
    try:
        Z_m=sorted(zip(np.vstack(m_HJD)-2.45e6,np.vstack(m_SiIII_fm)-np.mean(m_SiIII_fm)))
        m_HJD,m_fm=[Z_m[i][0] for i in range(len(Z_m))],[Z_m[i][1] for i in range(len(Z_m))]
        m_HJD,m_fm=np.vstack(m_HJD),np.vstack(m_fm)
    except: pass
    try:
        Z_f=sorted(zip(np.vstack(f_HJD)-2.45e6,np.vstack(f_SiIII_fm)-np.mean(f_SiIII_fm)))
        f_HJD,f_fm=[Z_f[i][0] for i in range(len(Z_f))],[Z_f[i][1] for i in range(len(Z_f))]
        f_HJDf_fm=np.vstack(f_HJD),np.vstack(f_fm)
    except: pass
    try:
        Z_s=sorted(zip(np.vstack(s_HJD)-2.45e6,np.vstack(s_SiIII_fm)-np.mean(s_SiIII_fm)))
        s_HJD,s_fm=[Z_s[i][0] for i in range(len(Z_s))],[Z_s[i][1] for i in range(len(Z_s))]
        s_HJD,s_fm=np.vstack(s_HJD),np.vstack(s_fm)
    except: pass

    STD_s_cv,STD_f_cv,STD_m_cv=[],[],[]
    #step = 40
    try:
        data_0=min(f_HJD[0],m_HJD[0],s_HJD[0])[0]
        data_1=max(f_HJD[-1],m_HJD[-1],s_HJD[-1])[0]
    except: 
        try:
            data_0=min(f_HJD[0],m_HJD[0])[0]
            data_1=max(f_HJD[-1],m_HJD[-1])[0]
        except: 
            try:
                data_0=min(f_HJD[0])[0]
                data_1=max(f_HJD[0])[0]
            except: pass
    #n_step=int((data_1-data_0)/step)+1
    for n in range(1,50):
        if len(f_HJD)!=0:
            STD_f_cv.append(np.nanstd(f_fm[int(0.02*n*len(f_fm)):int((0.02*n+0.1)*n*len(f_fm))]))
        if len(m_HJD)!=0:
            STD_m_cv.append(np.nanstd(m_fm[int(0.02*n*len(m_fm)):int((0.02*n+0.1)*len(m_fm))]))
        if len(s_HJD)!=0:
            STD_s_cv.append(np.nanstd(s_fm[int(0.02*n*len(s_fm)):int((0.02*n+0.1)*n*len(s_fm))])) 






    np.random.seed(5)
    STD_s_rn,STD_f_rn,STD_m_rn=[],[],[]
    #fig, axs=plt.subplots(sharex=True)
    for i in range(50000):
        # Selection of the data to be used:
        if len(s_HJD)!=0:
            num=[]
            rnd=np.random.randint(0,len(s_HJD),int(0.1*len(s_HJD)))
            for j in rnd:
                num.append(s_fm[j])
            STD_s_rn.append(np.std(num))
        if len(f_HJD)!=0:
            num=[]
            rnd=np.random.randint(0,len(f_HJD),int(0.1*len(f_HJD)))
            for j in rnd:
                num.append(f_fm[j])
            STD_f_rn.append(np.std(num))
        if len(m_HJD)!=0:
            num=[]
            rnd=np.random.randint(0,len(m_HJD),int(0.2*len(m_HJD)))   
            for j in rnd:
                num.append(m_fm[j])
            STD_m_rn.append(np.std(num))    

    pp_f,pp_m,pp_s=[],[],[]
    for i in range(25000):
        # Selection of the data to be used:
        if len(s_HJD)!=0:            
            num=[]
            rnd=np.random.randint(0,len(s_HJD),int(0.1*len(s_HJD)))
            for j in rnd:
                num.append(float(s_fm[j]))
            pp_s.append(max(num)-min(num))
        if len(f_HJD)!=0:            
            num=[]
            rnd=np.random.randint(0,len(f_HJD),int(0.1*len(f_HJD)))
            for j in rnd:
                num.append(float(f_fm[j]))
            pp_f.append(max(num)-min(num))
        if len(m_HJD)!=0:            
            num=[]
            rnd=np.random.randint(0,len(m_HJD),int(0.1*len(m_HJD)))
            for j in rnd:
                num.append(float(m_fm[j]))
            pp_m.append(max(num)-min(num)) 





    plt.close('all')
    try: plt.close(sp)
    except: pass
    sp,ax1=plt.subplots(2,2)
    if len(s_HJD)!=0:   
        ax1[0,0].plot(np.vstack(s_HJD),np.vstack(s_fm)-np.mean(s_fm),'.',label='Rescaled SONG')
        h_s=ax1[1,0].hist(pp_s,bins=30,density=True,alpha=0.3,label='SONG')
        ax1[1,0].plot([h_s[1][h_s[0].argmax()]-2*np.mean(error_s),h_s[1][h_s[0].argmax()]+2*np.mean(error_s)],[max(h_s[0]),max(h_s[0])],color='#1f77b4')
        ax1[1,0].plot([h_s[1][h_s[0].argmax()]-2*np.mean(error_s),h_s[1][h_s[0].argmax()]+2*np.mean(error_s)],[max(h_s[0]),max(h_s[0])],color='#1f77b4')
        rn_s=ax1[0,1].hist(STD_s_rn,bins=50,density=True,alpha=0.3,label='SONG')
        ax1[1,1].hist(STD_s_cv,bins=20,density=True,alpha=0.3,label='SONG')
        ax1[0,0].add_patch(Rectangle((data_0, rn_s[1][rn_s[0].argmax()]-np.std(rn_s[1])),data_1-data_0, 2*np.std(rn_s[1]),alpha=0.3,color='#1f77b4'))
        ax1[0,0].add_patch(Rectangle((data_0, -rn_s[1][rn_s[0].argmax()]-np.std(rn_s[1])),data_1-data_0, 2*np.std(rn_s[1]),alpha=0.3,color='#1f77b4'))
    if len(f_HJD)!=0:     
        ax1[0,0].plot(np.vstack(f_HJD),np.vstack(f_fm)-np.mean(f_fm),'.',label='Rescaled FIES')  
        h_f=ax1[1,0].hist(pp_f,bins=30,density=True,alpha=0.3,label='FIES')
        ax1[1,0].plot([h_f[1][h_f[0].argmax()]-2*np.mean(error_f),h_f[1][h_f[0].argmax()]+2*np.mean(error_f)],[max(h_f[0]),max(h_f[0])],color='#1f7f0e')
        ax1[1,0].plot([h_f[1][h_f[0].argmax()]-2*np.mean(error_f),h_f[1][h_f[0].argmax()]+2*np.mean(error_f)],[max(h_f[0]),max(h_f[0])],color='#1f7f0e')
        rn_f=ax1[0,1].hist(STD_f_rn,bins=50,density=True,alpha=0.3,label='FIES')
        ax1[1,1].hist(STD_f_cv,bins=20,density=True,alpha=0.3,label='FIES') 
        ax1[0,0].add_patch(Rectangle((data_0, rn_f[1][rn_f[0].argmax()]-np.std(rn_f[1])),data_1-data_0, 2*np.std(rn_f[1]),alpha=0.3,color='#1f7f0e'))
        ax1[0,0].add_patch(Rectangle((data_0, -rn_f[1][rn_f[0].argmax()]-np.std(rn_f[1])),data_1-data_0, 2*np.std(rn_f[1]),alpha=0.3,color='#1f7f0e'))
    if len(m_HJD)!=0: 
        ax1[0,0].plot(np.vstack(m_HJD),np.vstack(m_fm)-np.mean(m_fm),'.',label='Rescaled HERMES')    
        h_m=ax1[1,0].hist(pp_m,bins=30,density=True,alpha=0.3,label='HERMES')
        ax1[1,0].plot([h_m[1][h_m[0].argmax()]-2*np.mean(error_m),h_m[1][h_m[0].argmax()]+2*np.mean(error_m)],[max(h_m[0]),max(h_m[0])],color='#2ca02c')
        ax1[1,0].plot([h_m[1][h_m[0].argmax()]-2*np.mean(error_m),h_m[1][h_m[0].argmax()]+2*np.mean(error_m)],[max(h_m[0]),max(h_m[0])],color='#2ca02c')
        rn_m=ax1[0,1].hist(STD_m_rn,bins=50,density=True,alpha=0.3,label='HERMES')
        ax1[1,1].hist(STD_m_cv,bins=20,density=True,alpha=0.3,label='HERMES') 
        ax1[0,0].add_patch(Rectangle((data_0, rn_m[1][rn_m[0].argmax()]-np.std(rn_m[1])),data_1-data_0, 2*np.std(rn_m[1]),alpha=0.3,color='#2ca02c'))
        ax1[0,0].add_patch(Rectangle((data_0, -rn_m[1][rn_m[0].argmax()]-np.std(rn_m[1])),data_1-data_0, 2*np.std(rn_m[1]),alpha=0.3,color='#2ca02c'))   
    ax1[0,0].grid(),ax1[0,0].legend()
    ax1[0,0].set_xlabel('HJD [days]-2450000'),ax1[0,0].set_ylabel(r'$\langle v\rangle-\langle v_0\rangle$'),
    ax1[1,1].legend()
    ax1[1,1].set_xlabel(r'$\sigma_v$ [km/s]')
    ax1[1,1].set_title('Convolution')
    ax1[0,1].get_shared_x_axes().join(ax1[0,1], ax1[1,1])
    ax1[0,1].set_xticklabels([])
    ax1[0,1].legend()
    ax1[0,1].set_title('Random selection')
    ax1[1,0].legend()
    ax1[1,0].set_xlabel('Peak-to-peak [km/s]')
    plt.show()

def peak_peak(star,lin_fol):
    SM=reading_summary('../stars/'+star+'/'+lin_fol+'/00_'+star+'_'+lin_fol+'.txt')
    m_fm,f_fm,s_fm=[],[],[]
    #m_HJD,f_HJD,s_HJD=[],[],[]
    #error_m,error_f,error_s=[],[],[]
    try: 
        OF,LF,HJD,ZM=zero_reading(star,'H_alpha')
    except:
        OF,LF,HJD,ZM=zero_reading(star,'H_alpha_6562')
    for i in range(len(SM.first_moment)):
        if SM.Instrument[i]=="'MERCATOR'":
            m_fm.append(SM.first_moment[i])#,m_HJD.append(SM.HJD[i]),error_m.append(SM.error_fm[i])
        elif SM.Instrument[i]=="'FIES'":
            f_fm.append(SM.first_moment[i])#,f_HJD.append(SM.HJD[i]),error_f.append(SM.error_fm[i])
        elif SM.Instrument[i]=="'SONG'":
            s_fm.append(SM.first_moment[i])#,s_HJD.append(SM.HJD[i]),error_s.append(SM.error_fm[i])
        else:
            print('Unknown instrument')
    m_zm,f_zm,s_zm=[],[],[]
    for i in range(len(LF)):            
        if 'FIES' in LF[i]:   
            f_zm.append(ZM[i])
        if 'MERCATOR' in LF[i]:
            m_zm.append(ZM[i])    
        if 'SONG' in LF[i]:
            s_zm.append(ZM[i])
    if len(f_fm)!=0:
        f_fm=np.vstack(f_fm)
        f_zm=np.vstack(f_zm)        
    if len(m_fm)!=0:
        m_fm=np.vstack(f_fm)
        m_zm=np.vstack(f_zm)
    if len(s_fm)!=0:
        s_fm=np.vstack(f_fm)
        s_zm=np.vstack(f_zm)      
    pp_fm_f,pp_fm_m,pp_fm_s=[],[],[]      
    pp_zm_f,pp_zm_m,pp_zm_s=[],[],[]
    for i in range(25000):
        # Selection of the data to be used:
        if len(s_fm)>2:            
            num_fm=[]            
            num_zm=[]
            rnd=np.random.randint(0,len(s_fm),int(0.1*len(s_fm))+3)
            for j in rnd:
                num_fm.append(float(s_fm[j]))
                num_zm.append(float(s_zm[j]))
            pp_fm_s.append(max(num_fm)-min(num_fm))
            pp_zm_s.append(max(num_zm)-min(num_zm))            
        if len(f_fm)>2:            
            num_fm=[]            
            num_zm=[]
            rnd=np.random.randint(0,len(f_fm),int(0.1*len(f_fm))+3)
            for j in rnd:
                num_fm.append(float(f_fm[j]))
                num_zm.append(float(f_zm[j]))
            pp_fm_f.append(max(num_fm)-min(num_fm))
            pp_zm_f.append(max(num_zm)-min(num_zm))            
        if len(m_fm)>2:            
            num_fm=[]            
            num_zm=[]
            rnd=np.random.randint(0,len(m_fm),int(0.1*len(m_fm))+3)
            for j in rnd:
                num_fm.append(float(m_fm[j]))
                num_zm.append(float(m_zm[j]))
            pp_fm_m.append(max(num_fm)-min(num_fm))
            pp_zm_m.append(max(num_zm)-min(num_zm))
    hdu=open('../stars/'+star+'/peak_to_peak.txt','a')
    if len(m_fm)!=0:            
        mode_m_fm=(np.histogram(pp_fm_m)[1][np.argmax(np.histogram(pp_fm_m)[0])]+np.histogram(pp_fm_m)[1][np.argmax(np.histogram(pp_fm_m)[0])+1])/2
        mode_m_zm=(np.histogram(pp_zm_m)[1][np.argmax(np.histogram(pp_zm_m)[0])]+np.histogram(pp_zm_m)[1][np.argmax(np.histogram(pp_zm_m)[0])+1])/2
        std_m_fm=np.std(pp_fm_m)
        std_m_zm=np.std(pp_zm_m)
        info=(mode_m_fm,std_m_fm,mode_m_zm,std_m_zm)
        hdu.writelines(str(info))
    if len(s_fm)!=0:            
        mode_s_fm=(np.histogram(pp_fm_s)[1][np.argmax(np.histogram(pp_fm_s)[0])]+np.histogram(pp_fm_s)[1][np.argmax(np.histogram(pp_fm_s)[0])+1])/2
        mode_s_zm=(np.histogram(pp_zm_s)[1][np.argmax(np.histogram(pp_zm_s)[0])]+np.histogram(pp_zm_s)[1][np.argmax(np.histogram(pp_zm_s)[0])+1])/2
        std_s_fm=np.std(pp_fm_s)
        std_s_zm=np.std(pp_zm_s)
        info=(mode_s_fm,std_s_fm,mode_s_zm,std_s_zm)
        hdu.writelines(str(info))
    if len(f_fm)!=0:            
        mode_f_fm=(np.histogram(pp_fm_f)[1][np.argmax(np.histogram(pp_fm_f)[0])]+np.histogram(pp_fm_f)[1][np.argmax(np.histogram(pp_fm_f)[0])+1])/2
        mode_f_zm=(np.histogram(pp_zm_f)[1][np.argmax(np.histogram(pp_zm_f)[0])]+np.histogram(pp_zm_f)[1][np.argmax(np.histogram(pp_zm_f)[0])+1])/2
        std_f_fm=np.std(pp_fm_f)
        std_f_zm=np.std(pp_zm_f)
        info=(mode_f_fm,std_f_fm,mode_f_zm,std_f_zm)
        hdu.writelines(str(info))
    try: 
        HJD_ph,mag=photometry(star)
        num,pp=[],[]
        rnd=np.random.randint(0,len(mag),int(0.1*len(mag)))
        for j in rnd:
            num.append(float(mag[j]))
        pp.append(max(num)-min(num))
        mode=(np.histogram(pp)[1][np.argmax(np.histogram(pp)[0])]+np.histogram(pp)[1][np.argmax(np.histogram(pp)[0])+1])/2
        std=np.std(pp)
        info=(mode,std)
        hdu.writelines(str(info))
    except:
        print(star,' without photometry')
    hdu.close()
    
    #h_s=ax1[1,0].hist(pp_s,bins=30,density=True,alpha=0.3,label='SONG')