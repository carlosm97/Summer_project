#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 14:05:27 2021

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
import csv

# In[Statistics]
# It should be run again. 
STATISTIC=True
if STATISTIC==True:
    if os.path.exists('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/statistic.txt'):
        os.remove('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/statistic.txt')
    hdu=open('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/statistic.txt','a')
    
    # To know the spectral type: 
    #SB(star)['SP_TYPE'][0]
    
    hdu.writelines('Star, Photospheric sigma, Wind sigma,  Mean wind, Photometric sigma, len(FM), len(ZM)')
    hdu.writelines('\n')
    try:            
        hdu.writelines(str(statistical_study('HD2905','4567',4567.84,-125,150,-400,400,ph_name='Si III'))) # B1Ia
        hdu.writelines('\n')
    except: print('Error HD2905')
    try:                
        hdu.writelines(str(statistical_study('HD8065','6456',6456.383,-125,-35,-400,300,ph_name='Fe II'))) # B9Iab
        hdu.writelines('\n')
    except: print('Error HD8065')
    try:                
        hdu.writelines(str(statistical_study('HD21291','6347',6347.11,-100,100,-200,170,ph_name='Si II'))) # B9Ia
        hdu.writelines('\n')
    except: print('Error HD21291')
    try:                
        hdu.writelines(str(statistical_study('HD30614','5592',5592.252,-150,200,-600,600,ph_name='O III'))) # O9Ia
        hdu.writelines('\n')
    except: print('Error HD30614')
    try:                
        hdu.writelines(str(statistical_study('HD31327','4567',4567.84,-90,90,-390,400,ph_name='Si III'))) # B2.5Ib
        hdu.writelines('\n')
    except: print('Error HD31327')
    try:            
        hdu.writelines(str(statistical_study('HD34078','5592',5592.252,-30,120,-750,750,ph_name='O III'))) # O9.5V
        hdu.writelines('\n')
    except: print('Error HD34078')
    try:                
        hdu.writelines(str(statistical_study('HD34085','6347',6347.11,-80,110,-200,200,ph_name='Si II',TESS=False))) # B8Ia
        hdu.writelines('\n')    
    except: print('Error HD34085')
    try:                
        hdu.writelines(str(statistical_study('HD34656','5592',5592.252,-150,160,-400,400,ph_name='O III'))) # O7.5II(f)
        hdu.writelines('\n')
    except: print('Error HD34656')
    try:                
        hdu.writelines(str(statistical_study('HD36512','5592',5592.252,-35,55,-1000,1000,ph_name='O III'))) # O9.7V
        hdu.writelines('\n')
    except: print('Error HD36512')
    try:            
        hdu.writelines(str(statistical_study('HD36861','5592',5592.252,-75,175,-400,400,ph_name='O III'))) # O8III
        hdu.writelines('\n')
    except: print('Error HD36861')
    try:                
        hdu.writelines(str(statistical_study('HD37128','4567',4567.84,-150,200,-400,350,ph_name='Si III'))) # B0Ia
        hdu.writelines('\n')
    except: print('Error HD37128')
    try:                
        hdu.writelines(str(statistical_study('HD38771','4567',4567.84,-150,180,-400,500,ph_name='Si III'))) # B0.5Ia
        hdu.writelines('\n')
    except: print('Error HD38771')
    try:                
        hdu.writelines(str(statistical_study('HD47839','5592',5592.252,-100,130,-600,600,ph_name='O III'))) # O7V+B1.5/2V
        hdu.writelines('\n') # Binary system 
    except: print('Error HD47839')
    try:            
        hdu.writelines(str(statistical_study('HD91316','4567',4567.84,-70,180,-300,300,ph_name='Si III'))) # B1Iab
        hdu.writelines('\n')
    except: print('Error HD91316')
    try:                
        hdu.writelines(str(statistical_study('HD164353','4567',4567.84,-80,75,-600,600,ph_name='Si III',TESS=False))) # B5I
        hdu.writelines('\n')
    except: print('Error HD164353')
    try:                
        hdu.writelines(str(statistical_study('HD188001','5592',5592.252,-150,200,-1000,1000,ph_name='O III',TESS=False))) # O7.5Iabf
        hdu.writelines('\n')
    except: print('Error HD188001')
    try:            
        hdu.writelines(str(statistical_study('HD188209','5592',5592.252,-150,150,-350,230,ph_name='O III'))) # O9.5Iab
        hdu.writelines('\n')
    except: print('Error HD188209')
    try:            
        hdu.writelines(str(statistical_study('HD190603','4567',4567.84,-100,110,-300,400,ph_name='Si III'))) # B1.5Ia
        hdu.writelines('\n')
    except: print('Error HD190603')
    try:                
        hdu.writelines(str(statistical_study('HD198478','4567',4567.84,-100,130,-250,200,ph_name='Si III'))) # B4Ia
        hdu.writelines('\n')
    except: print('Error HD198478')
    try:                
        hdu.writelines(str(statistical_study('HD199478','4567',4567.84,-100,75,-400,300,ph_name='Si III'))) # B8Ia
        hdu.writelines('\n')
    except: print('Error HD199478')
    try:                
        hdu.writelines(str(statistical_study('HD206165','4567',4567.84,-110,110,-500,500,ph_name='Si III'))) # B2Ib
        hdu.writelines('\n')
    except: print('Error HD206165')
    try:                
        hdu.writelines(str(statistical_study('HD207198','5592',5592.252,-150,140,-400,400,ph_name='O III'))) # O8.5II
        hdu.writelines('\n')
    except: print('Error HD207198')
    try:                
        hdu.writelines(str(statistical_study('HD209975','4567',4567.84,-130,90,-500,500,ph_name='Si III'))) # O9Ib
        hdu.writelines('\n')
    except: print('Error HD209975')
    try:                
        hdu.writelines(str(statistical_study('HD213087','4567',4567.84,-150,150,-330,300,ph_name='Si III'))) # B0.5Ib
        hdu.writelines('\n')
    except: print('Error HD213087')
    try:                
        hdu.writelines(str(statistical_study('HD214680','5592',5592.252,-75,55,-750,750,ph_name='O III'))) # O9V
        hdu.writelines('\n')
    except: print('Error HD214680')
    
    hdu.close()
# In[final_figure]
star,folder,line,v1,v2,v_alpha_1,v_alpha_2='HD30614','5592',5592.252,-150,200,-600,600
plt.rcParams["figure.figsize"]=(8,6)
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
files=os.listdir(path)
common_fm=[]
#line=line
fm,HJD_fm=[],[]
if os.path.exists(path+'/uniform_fm.txt')==True:
    os.remove(path+'/uniform_fm.txt')
hdu=open(path+'/uniform_fm.txt','a')
X_int,Y_int=[],[]
if len(files)>55:
    criteria=50/len(files)
else: 
    criteria=1
M=[[],[]]
for file in files: 
    if file.startswith('HD') and file.endswith('.txt'):
        wv,v,f=reading_line_file(path+'/'+file)
        x_l,v_l,y_l=wv[v>v1],v[v>v1],f[v>v1]
        x_l,y_l,v_l=x_l[v_l<v2],y_l[v_l<v2],v_l[v_l<v2]
        llh=np.random.uniform(0,1)
        X_int.append(np.around(v,3)),Y_int.append(np.around(f,3))
        M[0].append(np.around(min(v),3)),M[1].append(np.around(max(v),3))
        FM=first_moment(x_l,y_l,line)/zero_moment(x_l,y_l,line)
        fm.append(FM)
        hjd=file.split('_')[0].replace(star,'')
        hjd=float(hjd.replace('p','.'))
        HJD_fm.append(hjd)
        if 'SONG' in file:
            ins='SONG'
        if 'MERCATOR' in file:
            ins='MERCATOR'
        if 'FIES' in file:
            ins='FIES'
        info=FM,hjd,ins
        info=str(info)
        hdu.writelines(info)
        hdu.writelines('\n')
X_int,Y_int=np.array(X_int),np.array(Y_int)
X_m,Y_m=[],[]
for j in range(len(X_int)):
    Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
    X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)    

if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
    os.mkdir('../stars/'+star+'/SUMMARY')

hdu=open(path+'/uniform_fm.txt','r')
FM_old,HJD_old_fm,INS_old_fm=[],[],[]
for j in hdu:
    dat=j.split()
    FM_old.append(float(dat[0][1:-1]))
    HJD_old_fm.append(float(dat[1][:-1]))
    INS_old_fm.append(dat[2][:-1])
FM_old=np.array(FM_old)-np.mean(FM_old)
HJD_old_fm=np.array(HJD_old_fm)-2.45e6
INS_old_fm=np.array(INS_old_fm)
FM_i,HJD_fm_i=[[],[],[]],[[],[],[]]
FM_r,HJD_fm_r=[[],[],[]],[[],[],[]]
hdu.close()
for fm in range(len(FM_old)):
    if 'MERCATOR' in INS_old_fm[fm]:
        ins=0
    elif 'FIES' in INS_old_fm[fm]:
        ins=1
    elif 'SONG' in INS_old_fm[fm]:
        ins=2
    if abs(FM_old[fm])<2.5*np.std(FM_old):
        FM_i[ins].append(FM_old[fm]),HJD_fm_i[ins].append(HJD_old_fm[fm])
    else:
        FM_r[ins].append(FM_old[fm]),HJD_fm_r[ins].append(HJD_old_fm[fm])
FM,HJD_fm=np.concatenate(FM_i),np.concatenate(HJD_fm_i)        
np.random.seed(5)

plt.close('all')
fig,axs=plt.subplots(3)
HJD_fm_1,FM_1=HJD_fm[HJD_fm<8100],FM[HJD_fm<8100]
HJD_fm_1,FM_1=HJD_fm_1[HJD_fm_1>8090],FM_1[HJD_fm_1>8090]
axs[0].plot(HJD_fm_1,FM_1,'.')
axs[0].grid()
axs[0].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[0].set_xlabel('HJD-2450000')

HJD_fm_2,FM_2=HJD_fm[HJD_fm<8841],FM[HJD_fm<8841]
HJD_fm_2,FM_2=HJD_fm_2[HJD_fm_2>8836],FM_2[HJD_fm_2>8836]
axs[2].plot(HJD_fm_2,FM_2,'.')
axs[2].grid()
axs[2].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[2].set_xlabel('HJD-2450000')

HJD_fm_3,FM_3=HJD_fm[HJD_fm<8110],FM[HJD_fm<8110]
HJD_fm_3,FM_3=HJD_fm_3[HJD_fm_3>8100],FM_3[HJD_fm_3>8100]
axs[1].plot(HJD_fm_3,FM_3,'.')
axs[1].grid()
axs[1].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[1].set_xlabel('HJD-2450000')
plt.show()
print(np.std(FM_1),np.std(FM_3),np.std(FM_2))
extra_8_ph=(np.std(FM_1)+np.std(FM_2)+np.std(FM_3))/3
star,folder,line,v1,v2='HD31327','4567',4567.84,-90,90
plt.rcParams["figure.figsize"]=(8,6)
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
files=os.listdir(path)
common_fm=[]
#line=line
fm,HJD_fm=[],[]
if os.path.exists(path+'/uniform_fm.txt')==True:
    os.remove(path+'/uniform_fm.txt')
hdu=open(path+'/uniform_fm.txt','a')
X_int,Y_int=[],[]
if len(files)>55:
    criteria=50/len(files)
else: 
    criteria=1
M=[[],[]]
for file in files: 
    if file.startswith('HD') and file.endswith('.txt'):
        wv,v,f=reading_line_file(path+'/'+file)
        x_l,v_l,y_l=wv[v>v1],v[v>v1],f[v>v1]
        x_l,y_l,v_l=x_l[v_l<v2],y_l[v_l<v2],v_l[v_l<v2]
        llh=np.random.uniform(0,1)
        X_int.append(np.around(v,3)),Y_int.append(np.around(f,3))
        M[0].append(np.around(min(v),3)),M[1].append(np.around(max(v),3))
        FM=first_moment(x_l,y_l,line)/zero_moment(x_l,y_l,line)
        fm.append(FM)
        hjd=file.split('_')[0].replace(star,'')
        hjd=float(hjd.replace('p','.'))
        HJD_fm.append(hjd)
        if 'SONG' in file:
            ins='SONG'
        if 'MERCATOR' in file:
            ins='MERCATOR'
        if 'FIES' in file:
            ins='FIES'
        info=FM,hjd,ins
        info=str(info)
        hdu.writelines(info)
        hdu.writelines('\n')
X_int,Y_int=np.array(X_int),np.array(Y_int)
X_m,Y_m=[],[]
for j in range(len(X_int)):
    Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
    X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)    

if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
    os.mkdir('../stars/'+star+'/SUMMARY')

hdu=open(path+'/uniform_fm.txt','r')
FM_old,HJD_old_fm,INS_old_fm=[],[],[]
for j in hdu:
    dat=j.split()
    FM_old.append(float(dat[0][1:-1]))
    HJD_old_fm.append(float(dat[1][:-1]))
    INS_old_fm.append(dat[2][:-1])
FM_old=np.array(FM_old)-np.mean(FM_old)
HJD_old_fm=np.array(HJD_old_fm)-2.45e6
INS_old_fm=np.array(INS_old_fm)
FM_i,HJD_fm_i=[[],[],[]],[[],[],[]]
FM_r,HJD_fm_r=[[],[],[]],[[],[],[]]
hdu.close()
for fm in range(len(FM_old)):
    if 'MERCATOR' in INS_old_fm[fm]:
        ins=0
    elif 'FIES' in INS_old_fm[fm]:
        ins=1
    elif 'SONG' in INS_old_fm[fm]:
        ins=2
    if abs(FM_old[fm])<2.5*np.std(FM_old):
        FM_i[ins].append(FM_old[fm]),HJD_fm_i[ins].append(HJD_old_fm[fm])
    else:
        FM_r[ins].append(FM_old[fm]),HJD_fm_r[ins].append(HJD_old_fm[fm])
FM,HJD_fm=np.concatenate(FM_i),np.concatenate(HJD_fm_i)        
np.random.seed(5)
plt.plot(HJD_fm,FM,'.')
plt.grid()

fig,axs=plt.subplots()
HJD_fm_1,FM_1=HJD_fm[HJD_fm<8922],FM[HJD_fm<8922]
HJD_fm_1,FM_1=HJD_fm_1[HJD_fm_1>8766],FM_1[HJD_fm_1>8766]
axs.plot(HJD_fm_1,FM_1,'.')
axs.grid()
axs.set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs.set_xlabel('HJD-2450000')
extra_18_ph=np.std(FM_1)






plt.rcParams["figure.figsize"] = (15,8)
plt.close('all')
hdu=open('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/statistic.txt','r')
s_ph,s_wnd,s_mg=[],[],[]
star=[]
ln_PH,ln_WND=[],[]
ns_st=[]
ST,LT=[],[]
fig,ax=plt.subplots(2,3)
for j in hdu: 
    try:
        j=j.replace('(','')
        j=j.replace(')','')
        dat=j.split(',')
        s_ph.append(float(dat[1])),s_wnd.append(abs(float(dat[2]))),s_mg.append(float(dat[4])) 
        ln_PH.append(float(dat[5])),ln_WND.append(float(dat[6]))
        star.append(dat[0][1:-1])
        ns_st.append(float(dat[0][3:-1]))
        spc=spc_code(SB(dat[0][1:-1])['SP_TYPE'][0])
        ST.append(spc[0]),LT.append(20*(7-spc[1]))
        #ax[0].plot(float(dat[3]),float(dat[1]),'.',markersize=10,label=dat[0])
        #ax[1].plot(float(dat[3]),float(dat[2]),'.',markersize=10,label=dat[0])
        #ax[2].plot(float(dat[1]),float(dat[2]),'.',markersize=10,label=dat[0])
        if float(dat[2])>3:
            print(dat[0][1:-1])
    except: 
        pass
new_star_order=['HD47839','HD214680','HD34078','HD36512','HD36861','HD188001',\
                'HD34656','HD207198','HD30614','HD209975','HD188209','HD37128',\
                'HD38771','HD213087','HD2905','HD91316','HD190603','HD206165',\
                'HD31327','HD198478','HD164353','HD199478','HD34085','HD21291','HD8065']    
order = {a: i for i, a in enumerate(new_star_order)}
odr=sorted(zip(star,s_ph,s_wnd,s_mg,ln_PH,ln_WND),key=lambda x: order[x[0]]) #sorted(zip(star,s_ph,s_wnd,s_mg,ln_PH,ln_WND))
star,s_ph=[odr[i][0] for i in range(0,len(star))],[odr[i][1] for i in range(0,len(star))]
s_wnd,s_mg=[odr[i][2] for i in range(0,len(star))],[odr[i][3] for i in range(0,len(star))]
ln_PH,ln_WND=[odr[i][4] for i in range(0,len(star))],[odr[i][5] for i in range(0,len(star))]

s_ph,s_wnd,s_mg=np.array(s_ph),np.array(s_wnd),np.array(s_mg)
hdu.close()
hdu=open('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/stars.txt','r')
STR,T_eff,L=[],[],[]
for j in hdu: 
    dat=j.split()
    if dat[0] in star:
        T_eff.append(4+np.log10(float(dat[1]))),L.append(5.39-float(dat[2])),STR.append(dat[0])
T_eff.append(4.32),L.append(4.33),STR.append('HD190603')
hdu.close()
logT,logL=[],[]
bc=[]
for j in star:
    try: 
        logT.append(np.array(T_eff)[np.array(STR)==j][0]),logL.append(np.array(L)[np.array(STR)==j][0]),
    except: logT.append(np.nan),logL.append(np.nan)
logL,logT=np.array(logL),np.array(logT)
L_g=50*(np.array(logL)-min(logL))+5    
alp=np.log10(ln_PH)/max(np.log10(ln_PH))
cm = plt.cm.get_cmap('brg_r')
numeration=range(1,len(s_ph)+1)

ax[0,0].plot([s_mg[8],s_mg[8]],[s_ph[8],extra_8_ph],'indigo',alpha=alp[8])
im=ax[0,0].scatter(s_mg,s_ph,c=logT,s=L_g,cmap=cm,alpha=alp)
ax[0,0].scatter(s_mg[8],extra_8_ph,c='indigo',s=L_g[8],marker="s",alpha=alp[8])
ax[0,0].plot([s_mg[18],s_mg[18]],[s_ph[18],extra_18_ph],'red',alpha=alp[18])
ax[0,0].scatter(s_mg[18],extra_18_ph,c='red',s=L_g[18],marker="s",alpha=alp[18])

ax[0,0].scatter(s_mg[np.isnan(logT)],s_ph[np.isnan(logT)],c='k')
alp=np.log10(ln_WND)/max(np.log10(ln_WND))
ax[0,2].plot([s_ph[8],extra_8_ph],[s_wnd[8],s_wnd[8]],'indigo')
ax[0,1].scatter(s_mg,s_wnd,c=logT,s=L_g,cmap=cm,alpha=alp)
ax[0,1].scatter(s_mg[np.isnan(logT)],s_wnd[np.isnan(logT)],c='k')
im2=ax[0,2].scatter(s_ph,s_wnd,c=logT,s=L_g,cmap=cm)
ax[0,2].scatter(extra_8_ph,s_wnd[8],s=L_g[8],marker="s",c='indigo')
ax[0,2].scatter(extra_18_ph,s_wnd[18],s=L_g[18],marker="s",c='red')
ax[0,2].plot([s_ph[18],extra_18_ph],[s_wnd[18],s_wnd[18]],'red')

ax[0,2].scatter(s_ph[np.isnan(logT)],s_wnd[np.isnan(logT)],c='k')
cb_ax = fig.add_axes([0.91, 0.11, 0.02, 0.77])
cbar=fig.colorbar(im2,cax=cb_ax)
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'$\log(T_{eff})$', rotation=270,labelpad=15,fontsize=14)
ax[0,0].get_shared_x_axes().join(ax[0,0], ax[0,1])
ax[0,1].get_shared_y_axes().join(ax[0,1], ax[0,2])
ax[0,0].grid()#,ax[0].legend()
ax[0,0].set_xlabel(r'TESS, $\sigma_{Tp}$ [mmag]',fontsize=14),ax[0,0].set_ylabel(r'Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]',fontsize=14)
ax[0,1].grid()#,ax[1].legend()
ax[0,1].set_xlabel(r'TESS, $\sigma_{Tp}$ [mmag]',fontsize=14),ax[0,1].set_ylabel(r'$H\alpha$, $\sigma_{\langle v^0\rangle}$ [km/s]', labelpad=-5,fontsize=14)
ax[0,2].grid()#,ax[0].legend()
#ax[0].legend(framealpha=0.1,loc='upper right')
ax[0,2].set_xlabel(r'Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]',fontsize=14),ax[0,2].set_ylabel(r'$H\alpha$, $\sigma_{\langle v^0\rangle}$ [km/s]', labelpad=-5,fontsize=14)
plt.show()
#ax[0,2].set_yscale("log"),ax[0,1].set_yscale("log")
for i, txt in enumerate(numeration):
    ax[0,0].annotate(txt, (s_mg[i],s_ph[i]+0.2))
    ax[0,1].annotate(txt, (s_mg[i]+0.5,s_wnd[i]))
    ax[0,2].annotate(txt, (s_ph[i]+0.1,s_wnd[i]))
    

s_ph[8]=extra_8_ph
s_ph[18]=extra_18_ph
ax[1,0].scatter(logT,logL,c=logT,s=5*s_mg,cmap=cm)
ax[1,1].scatter(logT,logL,c=logT,s=5*s_wnd,cmap=cm)
ax[1,2].scatter(logT,logL,c=logT,s=10*s_ph,cmap=cm)



path='/home/charlie/Desktop/MSc/Primero/Primer cuatrimetre/Estructuras y evolucion estelar/Eassy/Modelos_Sergio'
modelos = os.listdir(path+'/modelos/')
modelos.sort()
# If you don't remember the name of any object you want to plot, uncomment the
# next lines to check an example. 
ZAMS=[[],[]]
ax[1,0].invert_xaxis(),ax[1,1].invert_xaxis(),ax[1,2].invert_xaxis()
M0,LOG_L=[],[]
logTc,logrhoc=[],[]
for mod in modelos:
    if mod[-5]=='0': #or mod[-5]=='4':
        if mod[1:4]=='009' or  mod[1:4]=='015' or mod[1:4]=='025'\
        or mod[1:4]=='040' or mod[1:4]=='085' or mod[1:4]=='120':
            Fich=lectura(path+'/modelos/'+mod)
            Log_L,Log_Teff=select(Fich,'lg(L)'),select(Fich,'lg(Teff)')
            mass=select(Fich,'mass')
            log_L_g=Log_L-np.log10(mass)
            ZAMS[0].append(Log_Teff[0])
            ZAMS[1].append(log_L_g[0])
            ax[1,0].plot(Log_Teff[:200],log_L_g[:200],c='k')
            ax[1,1].plot(Log_Teff[:200],log_L_g[:200],c='k')
            ax[1,2].plot(Log_Teff[:200],log_L_g[:200],c='k')
            if float(mod[1:4])<80 and float(mod[1:4])>9:
                ax[1,0].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
                ax[1,1].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
                ax[1,2].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
            if float(mod[1:4])==85:
                ax[1,0].text(4.63,log_L_g[60]-0.0,str(mod[2:4])+r'$M_{\odot}$')
                ax[1,1].text(4.63,log_L_g[60]-0.0,str(mod[2:4])+r'$M_{\odot}$')
                ax[1,2].text(4.63,log_L_g[60]-0.,str(mod[2:4])+r'$M_{\odot}$')
            if float(mod[1:4])==120:
                ax[1,0].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
                ax[1,1].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
                ax[1,2].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
            #plt.xlabel(r'log($T_{eff}$)')
            #plt.ylabel('log(L)')
            #plt.legend()
            #plt.title('H-R with initial rotation velocity 0')
ZAMS[0].sort(),ZAMS[1].sort()
ax[1,0].plot(ZAMS[0],ZAMS[1],'k')
ax[1,0].grid()
ax[1,0].set_xlabel(r'log($T_{eff}$[K])',fontsize=14),ax[1,0].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)',fontsize=14)
ax[1,0].set_ylim(min(logL)-0.1,max(logL)+0.1)
ax[1,1].plot(ZAMS[0],ZAMS[1],'k')
ax[1,1].grid()
ax[1,1].set_xlabel(r'log($T_{eff}$[K])',fontsize=14),ax[1,1].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)',fontsize=14)
ax[1,1].set_ylim(min(logL)-0.1,max(logL)+0.1)
ax[1,2].plot(ZAMS[0],ZAMS[1],'k')
ax[1,2].grid()
ax[1,2].set_xlabel(r'log($T_{eff}$[K])',fontsize=14),ax[1,2].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)',fontsize=14)
ax[1,2].set_ylim(min(logL)-0.1,max(logL)+0.1)
if os.path.exists('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt')==True:
    os.remove('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt')
hdu=open('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt','w')
for i, txt in enumerate(numeration):
    if txt!=6 and txt!=21 and txt!=23 :
        ax[1,0].annotate(txt, (logT[i]-0.02,logL[i]))
    ax[1,1].annotate(txt, (logT[i]-0.02,logL[i]))
    ax[1,2].annotate(txt, (logT[i]-0.02,logL[i]))
    hdu.writelines(star[i])
    hdu.writelines(' ')
    hdu.writelines(str(txt))
    hdu.writelines('\n')
hdu.close()
ax[1,0].plot([],[],' ',label=r'size: TESS, $\sigma_{Tp}$ [mmag]')
ax[1,1].plot([],[],' ',label=r'size: $H\alpha$, $\sigma_{\langle v^0\rangle}$ [km/s]')
ax[1,2].plot([],[],' ',label=r'size: Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]')
ax[1,0].legend(framealpha=0.1,loc=4),ax[1,1].legend(framealpha=0.1,loc=4),ax[1,2].legend(framealpha=0.1,loc=4)
ax[1,0].set_xlim(4.73,3.95),ax[1,1].set_xlim(4.73,3.95),ax[1,2].set_xlim(4.73,3.95)
ax[0,0].tick_params(axis='y', labelsize=14),ax[0,0].tick_params(axis='x', labelsize=14) 
ax[0,1].tick_params(axis='y', labelsize=14),ax[0,1].tick_params(axis='x', labelsize=14)
ax[0,2].tick_params(axis='y', labelsize=14),ax[0,2].tick_params(axis='x', labelsize=14)
ax[1,0].tick_params(axis='y', labelsize=14),ax[1,0].tick_params(axis='x', labelsize=14)
ax[1,1].tick_params(axis='y', labelsize=14),ax[1,1].tick_params(axis='x', labelsize=14)
ax[1,2].tick_params(axis='y', labelsize=14),ax[1,2].tick_params(axis='x', labelsize=14)

plt.savefig('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/SUMMARY/HR_sigma.png')
# Mark with circles the studied stars:

num='HD188209','HD30614','HD2905','HD198478','HD190603'
for i in range(len(star)):
    if star[i] in num:
        ax[0,0].scatter(s_mg[i],s_ph[i],facecolors='none',edgecolors='k',s=100)#L_g[i-1]+5)
        ax[0,1].scatter(s_mg[i],s_wnd[i],facecolors='none',edgecolors='k',s=100)
        ax[0,2].scatter(s_ph[i],s_wnd[i],facecolors='none',edgecolors='k',s=100)
        ax[1,0].scatter(logT[i],logL[i],facecolors='none',edgecolors='k',s=150)
        ax[1,1].scatter(logT[i],logL[i],facecolors='none',edgecolors='k',s=150)
        ax[1,2].scatter(logT[i],logL[i],facecolors='none',edgecolors='k',s=150)
        
 
        
        
        
plt.savefig('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/SUMMARY/HR_sigma_marked.png')
        
plt.figure()
plt.scatter(s_ph,s_wnd,c=logT,s=L_g,cmap=cm,alpha=alp)
plt.grid()
# In[GIF]
H_alpha=6562.80 
GIF=False
if GIF==True:
    plt.close('all')
    star='HD198478'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,None)
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2457570,2457630,ll1=8710,ll2=8765)
    plt.show()
    
    
    plt.close('all')
    star='HD37128'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,'SONG')
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2457375,2457385,ll1=9170,ll2=9205,pbt=0.2)
    plt.show()
    plt.close('all')
    
    plt.close('all') 
    star='HD188209'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha',H_alpha,'SONG')
    gif(star,'H_alpha',HJD,F,V,Fl_Mean,2457590,2457620,ll1=8680,ll2=8740)
    plt.show()
    plt.close('all')
    
    plt.close('all') 
    star='HD30614'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,'FIES')
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2458836,2458841)
    plt.show()
    plt.close('all')
    
    plt.close('all') 
    star='HD2905'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,'SONG')
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2457947,2458077.5,ll1=8764,ll2=8815.5)
    plt.show()
    plt.close('all')
    
    '''
    plt.close('all')
    star='HD8065'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,None)
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2458768,2458922,ll1=8791,ll2=8841)
    plt.show()
    plt.close('all')
    '''
    plt.close('all') 
    star='HD190603'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,None)
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2458683,2458711)
    plt.show()
    plt.close('all')
    
    plt.close('all') 
    star='HD31327'
    WV,V,F,HJD,Fl_Mean=mean_line(star,'H_alpha_6562',H_alpha,None)
    gif(star,'H_alpha_6562',HJD,F,V,Fl_Mean,2458815,2458843)
    plt.show()
    plt.close('all')
# In[Image_copy]
# Image copy to different folders:
import shutil
stars=os.listdir('../stars')
shutil. rmtree('../stars/SUMMARY/H_alpha')
shutil. rmtree('../stars/SUMMARY/Photosphere')
shutil. rmtree('../stars/SUMMARY/Statistical_summary')
os.mkdir('../stars/SUMMARY/H_alpha')
os.mkdir('../stars/SUMMARY/Photosphere')
os.mkdir('../stars/SUMMARY/Statistical_summary')
for fld in stars: 
    if fld.startswith('HD'):
        if os.path.exists('../stars/'+fld+'/SUMMARY'):
            try:
                shutil.copy('../stars/'+fld+'/SUMMARY/'+fld+'_wind.png','../stars/SUMMARY/H_alpha/'+fld+'_wind.png')
                shutil.copy('../stars/'+fld+'/SUMMARY/'+fld+'_photosphere.png','../stars/SUMMARY/Photosphere/'+fld+'_photosphere.png')
                shutil.copy('../stars/'+fld+'/SUMMARY/'+fld+'_all.png','../stars/SUMMARY/Statistical_summary/'+fld+'_all.png')
            except: 
                print(fld)
'''
#%%
plt.close('all') 
star='HD31327'
fig,ax=plt.subplots()
WV,V,F,HJD,Fl_Mean=mean_line(star,'4567',4567.84,'FIES')
for i in range(len(WV)):
    ax.plot(V[i],F[i],'b',alpha=0.1)
ax.plot(V[0],Fl_Mean,'b')
WV,V,F,HJD,Fl_Mean=mean_line(star,'4567',4567.84,'MERCATOR')
for i in range(len(WV)):
    ax.plot(V[i],F[i],'g',alpha=0.1)
ax.plot(V[0],Fl_Mean,'g')
WV,V,F,HJD,Fl_Mean=mean_line(star,'4567',4567.84,'SONG')
for i in range(len(WV)):
    if np.random.rand()<0.1:
        ax.plot(V[i],F[i],'r',alpha=0.05)
ax.plot(V[0],Fl_Mean,'r')
plt.grid()

prueba_linea(star,'FIES',4567.84,4562.84,4572.84,x='vel')
prueba_linea(star,'SONG',4567.84,4562.84,4572.84,x='vel')

prueba_linea(star,'FIES',H_alpha,H_alpha-3,H_alpha+3,x='vel')
prueba_linea(star,'SONG',H_alpha,H_alpha-3,H_alpha+3,x='vel')
'''
# In[Window_study_30614]
star,folder,line,v1,v2,v_alpha_1,v_alpha_2='HD30614','5592',5592.252,-150,200,-600,600
plt.rcParams["figure.figsize"]=(8,6)
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
files=os.listdir(path)
common_fm=[]
#line=line
fm,HJD_fm=[],[]
if os.path.exists(path+'/uniform_fm.txt')==True:
    os.remove(path+'/uniform_fm.txt')
hdu=open(path+'/uniform_fm.txt','a')
X_int,Y_int=[],[]
if len(files)>55:
    criteria=50/len(files)
else: 
    criteria=1
M=[[],[]]
for file in files: 
    if file.startswith('HD') and file.endswith('.txt'):
        wv,v,f=reading_line_file(path+'/'+file)
        x_l,v_l,y_l=wv[v>v1],v[v>v1],f[v>v1]
        x_l,y_l,v_l=x_l[v_l<v2],y_l[v_l<v2],v_l[v_l<v2]
        llh=np.random.uniform(0,1)
        if llh<criteria:
            axfm.plot(v,f,alpha=0.5)
        X_int.append(np.around(v,3)),Y_int.append(np.around(f,3))
        M[0].append(np.around(min(v),3)),M[1].append(np.around(max(v),3))
        FM=first_moment(x_l,y_l,line)/zero_moment(x_l,y_l,line)
        fm.append(FM)
        hjd=file.split('_')[0].replace(star,'')
        hjd=float(hjd.replace('p','.'))
        HJD_fm.append(hjd)
        if 'SONG' in file:
            ins='SONG'
        if 'MERCATOR' in file:
            ins='MERCATOR'
        if 'FIES' in file:
            ins='FIES'
        info=FM,hjd,ins
        info=str(info)
        hdu.writelines(info)
        hdu.writelines('\n')
X_int,Y_int=np.array(X_int),np.array(Y_int)
X_m,Y_m=[],[]
for j in range(len(X_int)):
    Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
    X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)    

if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
    os.mkdir('../stars/'+star+'/SUMMARY')

hdu=open(path+'/uniform_fm.txt','r')
FM_old,HJD_old_fm,INS_old_fm=[],[],[]
for j in hdu:
    dat=j.split()
    FM_old.append(float(dat[0][1:-1]))
    HJD_old_fm.append(float(dat[1][:-1]))
    INS_old_fm.append(dat[2][:-1])
FM_old=np.array(FM_old)-np.mean(FM_old)
HJD_old_fm=np.array(HJD_old_fm)-2.45e6
INS_old_fm=np.array(INS_old_fm)
FM_i,HJD_fm_i=[[],[],[]],[[],[],[]]
FM_r,HJD_fm_r=[[],[],[]],[[],[],[]]
hdu.close()
for fm in range(len(FM_old)):
    if 'MERCATOR' in INS_old_fm[fm]:
        ins=0
    elif 'FIES' in INS_old_fm[fm]:
        ins=1
    elif 'SONG' in INS_old_fm[fm]:
        ins=2
    if abs(FM_old[fm])<2.5*np.std(FM_old):
        FM_i[ins].append(FM_old[fm]),HJD_fm_i[ins].append(HJD_old_fm[fm])
    else:
        FM_r[ins].append(FM_old[fm]),HJD_fm_r[ins].append(HJD_old_fm[fm])
FM,HJD_fm=np.concatenate(FM_i),np.concatenate(HJD_fm_i)        
np.random.seed(5)

plt.close('all')
fig,axs=plt.subplots(3)
HJD_fm_1,FM_1=HJD_fm[HJD_fm<8100],FM[HJD_fm<8100]
HJD_fm_1,FM_1=HJD_fm_1[HJD_fm_1>8090],FM_1[HJD_fm_1>8090]
axs[0].plot(HJD_fm_1,FM_1,'.')
axs[0].grid()
axs[0].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[0].set_xlabel('HJD-2450000')

HJD_fm_2,FM_2=HJD_fm[HJD_fm<8841],FM[HJD_fm<8841]
HJD_fm_2,FM_2=HJD_fm_2[HJD_fm_2>8836],FM_2[HJD_fm_2>8836]
axs[2].plot(HJD_fm_2,FM_2,'.')
axs[2].grid()
axs[2].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[2].set_xlabel('HJD-2450000')

HJD_fm_3,FM_3=HJD_fm[HJD_fm<8110],FM[HJD_fm<8110]
HJD_fm_3,FM_3=HJD_fm_3[HJD_fm_3>8100],FM_3[HJD_fm_3>8100]
axs[1].plot(HJD_fm_3,FM_3,'.')
axs[1].grid()
axs[1].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[1].set_xlabel('HJD-2450000')
plt.show()
print(np.std(FM_1),np.std(FM_3),np.std(FM_2))
'''
star,folder,line,v1,v2,v_alpha_1,v_alpha_2='HD190603','4567',4567.84,-100,110,-300,400
plt.rcParams["figure.figsize"]=(8,6)
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
files=os.listdir(path)
common_fm=[]
#line=line
fm,HJD_fm=[],[]
if os.path.exists(path+'/uniform_fm.txt')==True:
    os.remove(path+'/uniform_fm.txt')
hdu=open(path+'/uniform_fm.txt','a')
X_int,Y_int=[],[]
if len(files)>55:
    criteria=50/len(files)
else: 
    criteria=1
M=[[],[]]
for file in files: 
    if file.startswith('HD') and file.endswith('.txt'):
        wv,v,f=reading_line_file(path+'/'+file)
        x_l,v_l,y_l=wv[v>v1],v[v>v1],f[v>v1]
        x_l,y_l,v_l=x_l[v_l<v2],y_l[v_l<v2],v_l[v_l<v2]
        llh=np.random.uniform(0,1)
        if llh<criteria:
            axfm.plot(v,f,alpha=0.5)
        X_int.append(np.around(v,3)),Y_int.append(np.around(f,3))
        M[0].append(np.around(min(v),3)),M[1].append(np.around(max(v),3))
        FM=first_moment(x_l,y_l,line)/zero_moment(x_l,y_l,line)
        fm.append(FM)
        hjd=file.split('_')[0].replace(star,'')
        hjd=float(hjd.replace('p','.'))
        HJD_fm.append(hjd)
        if 'SONG' in file:
            ins='SONG'
        if 'MERCATOR' in file:
            ins='MERCATOR'
        if 'FIES' in file:
            ins='FIES'
        info=FM,hjd,ins
        info=str(info)
        hdu.writelines(info)
        hdu.writelines('\n')
X_int,Y_int=np.array(X_int),np.array(Y_int)
X_m,Y_m=[],[]
for j in range(len(X_int)):
    Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
    X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)    

if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
    os.mkdir('../stars/'+star+'/SUMMARY')

hdu=open(path+'/uniform_fm.txt','r')
FM_old,HJD_old_fm,INS_old_fm=[],[],[]
for j in hdu:
    dat=j.split()
    FM_old.append(float(dat[0][1:-1]))
    HJD_old_fm.append(float(dat[1][:-1]))
    INS_old_fm.append(dat[2][:-1])
FM_old=np.array(FM_old)-np.mean(FM_old)
HJD_old_fm=np.array(HJD_old_fm)-2.45e6
INS_old_fm=np.array(INS_old_fm)
FM_i,HJD_fm_i=[[],[],[]],[[],[],[]]
FM_r,HJD_fm_r=[[],[],[]],[[],[],[]]
hdu.close()
for fm in range(len(FM_old)):
    if 'MERCATOR' in INS_old_fm[fm]:
        ins=0
    elif 'FIES' in INS_old_fm[fm]:
        ins=1
    elif 'SONG' in INS_old_fm[fm]:
        ins=2
    if abs(FM_old[fm])<2.5*np.std(FM_old):
        FM_i[ins].append(FM_old[fm]),HJD_fm_i[ins].append(HJD_old_fm[fm])
    else:
        FM_r[ins].append(FM_old[fm]),HJD_fm_r[ins].append(HJD_old_fm[fm])
FM,HJD_fm=np.concatenate(FM_i),np.concatenate(HJD_fm_i)        
np.random.seed(5)

plt.close('all')
fig,axs=plt.subplots(2)
HJD_fm_1,FM_1=HJD_fm[HJD_fm<7750],FM[HJD_fm<7750]
HJD_fm_1,FM_1=HJD_fm_1[HJD_fm_1>7600],FM_1[HJD_fm_1>7600]
axs[0].plot(HJD_fm_1,FM_1,'.')
axs[0].grid()
axs[0].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[0].set_xlabel('HJD-2450000')
HJD_fm_2,FM_2=HJD_fm[HJD_fm<8950],FM[HJD_fm<8950]
HJD_fm_2,FM_2=HJD_fm_2[HJD_fm_2>8650],FM_2[HJD_fm_2>8650]
axs[1].plot(HJD_fm_2,FM_2,'.')
axs[1].grid()
axs[1].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[1].set_xlabel('HJD-2450000')
plt.show()
print(np.std(FM_1),np.std(FM_2))
'''
star=['HD47839','HD214680','HD34078','HD36512','HD36861','HD188001',\
                'HD34656','HD207198','HD30614','HD209975','HD188209','HD37128',\
                'HD38771','HD213087','HD2905','HD91316','HD190603','HD206165',\
                'HD31327','HD198478','HD164353','HD199478','HD34085','HD21291','HD8065'] 
s_ph[8]=(np.std(FM_1)+np.std(FM_2)+np.std(FM_3))/3

fig,ax=plt.subplots(2,3)
numeration=range(1,len(s_ph)+1)

im=ax[0,0].scatter(s_mg,s_ph,c=logT,s=L_g,cmap=cm,alpha=alp)
ax[0,0].scatter(s_mg[np.isnan(logT)],s_ph[np.isnan(logT)],c='k')
alp=np.log10(ln_WND)/max(np.log10(ln_WND))
ax[0,1].scatter(s_mg,s_wnd,c=logT,s=L_g,cmap=cm,alpha=alp)
ax[0,1].scatter(s_mg[np.isnan(logT)],s_wnd[np.isnan(logT)],c='k')
im2=ax[0,2].scatter(s_ph,s_wnd,c=logT,s=L_g,cmap=cm)
ax[0,2].scatter(s_ph[np.isnan(logT)],s_wnd[np.isnan(logT)],c='k')
cb_ax = fig.add_axes([0.91, 0.11, 0.02, 0.77])
cbar=fig.colorbar(im2,cax=cb_ax)
cbar.set_label(r'$\log(T_{eff})$', rotation=270,labelpad=15,fontasize=)
ax[0,0].get_shared_x_axes().join(ax[0,0], ax[0,1])
ax[0,1].get_shared_y_axes().join(ax[0,1], ax[0,2])
ax[0,0].grid()#,ax[0].legend()
ax[0,0].set_xlabel(r'TESS, $\sigma_{Tp}$ [mmag]'),ax[0,0].set_ylabel(r'Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]')
ax[0,1].grid()#,ax[1].legend()
ax[0,1].set_xlabel(r'TESS, $\sigma_{Tp}$ [mmag]'),ax[0,1].set_ylabel(r'$H\alpha$, $\sigma_{\langle v^0\rangle}$', labelpad=-5)
ax[0,2].grid()#,ax[0].legend()
#ax[0].legend(framealpha=0.1,loc='upper right')
ax[0,2].set_xlabel(r'Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]'),ax[0,2].set_ylabel(r'$H\alpha$, $\sigma_{\langle v^0\rangle}$', labelpad=-5)
plt.show()
ax[0,2].set_yscale("log"),ax[0,1].set_yscale("log")
for i, txt in enumerate(numeration):
    ax[0,0].annotate(txt, (s_mg[i],s_ph[i]+0.2))
    ax[0,1].annotate(txt, (s_mg[i]+0.5,s_wnd[i]))
    ax[0,2].annotate(txt, (s_ph[i]+0.1,s_wnd[i]))
    


ax[1,0].scatter(logT,logL,c=logT,s=5*s_mg,cmap=cm)
ax[1,1].scatter(logT,logL,c=logT,s=50*s_wnd,cmap=cm)
ax[1,2].scatter(logT,logL,c=logT,s=10*s_ph,cmap=cm)



path='/home/charlie/Desktop/MSc/Primero/Primer cuatrimetre/Estructuras y evolucion estelar/Eassy/Modelos_Sergio'
modelos = os.listdir(path+'/modelos/')
modelos.sort()
# If you don't remember the name of any object you want to plot, uncomment the
# next lines to check an example. 
ZAMS=[[],[]]
ax[1,0].invert_xaxis(),ax[1,1].invert_xaxis(),ax[1,2].invert_xaxis()
M0,LOG_L=[],[]
logTc,logrhoc=[],[]
for mod in modelos:
    if mod[-5]=='0': #or mod[-5]=='4':
        if mod[1:4]=='009' or  mod[1:4]=='015' or mod[1:4]=='025'\
        or mod[1:4]=='040' or mod[1:4]=='085' or mod[1:4]=='120':
            Fich=lectura(path+'/modelos/'+mod)
            Log_L,Log_Teff=select(Fich,'lg(L)'),select(Fich,'lg(Teff)')
            mass=select(Fich,'mass')
            log_L_g=Log_L-np.log10(mass)
            ZAMS[0].append(Log_Teff[0])
            ZAMS[1].append(log_L_g[0])
            ax[1,0].plot(Log_Teff[:200],log_L_g[:200],c='k')
            ax[1,1].plot(Log_Teff[:200],log_L_g[:200],c='k')
            ax[1,2].plot(Log_Teff[:200],log_L_g[:200],c='k')
            if float(mod[1:4])<80 and float(mod[1:4])>9:
                ax[1,0].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
                ax[1,1].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
                ax[1,2].text(4.15,log_L_g[190],str(mod[2:4])+r'$M_{\odot}$')
            if float(mod[1:4])==85:
                ax[1,0].text(4.63,log_L_g[60]-0.0,str(mod[2:4])+r'$M_{\odot}$')
                ax[1,1].text(4.63,log_L_g[60]-0.0,str(mod[2:4])+r'$M_{\odot}$')
                ax[1,2].text(4.63,log_L_g[60]-0.,str(mod[2:4])+r'$M_{\odot}$')
            if float(mod[1:4])==120:
                ax[1,0].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
                ax[1,1].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
                ax[1,2].text(4.7,log_L_g[60]-0.0,str(mod[1:4])+r'$M_{\odot}$')
            #plt.xlabel(r'log($T_{eff}$)')
            #plt.ylabel('log(L)')
            #plt.legend()
            #plt.title('H-R with initial rotation velocity 0')
ZAMS[0].sort(),ZAMS[1].sort()
ax[1,0].plot(ZAMS[0],ZAMS[1],'k')
ax[1,0].grid()
ax[1,0].set_xlabel(r'log($T_{eff}$[K])'),ax[1,0].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)')
ax[1,0].set_ylim(min(logL)-0.1,max(logL)+0.1)
ax[1,1].plot(ZAMS[0],ZAMS[1],'k')
ax[1,1].grid()
ax[1,1].set_xlabel(r'log($T_{eff}$[K])'),ax[1,1].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)')
ax[1,1].set_ylim(min(logL)-0.1,max(logL)+0.1)
ax[1,2].plot(ZAMS[0],ZAMS[1],'k')
ax[1,2].grid()
ax[1,2].set_xlabel(r'log($T_{eff}$[K])'),ax[1,2].set_ylabel(r'log($\mathcal{L/L}_{\odot}$)')
ax[1,2].set_ylim(min(logL)-0.1,max(logL)+0.1)
if os.path.exists('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt')==True:
    os.remove('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt')
hdu=open('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_study/star_correspondence.txt','w')
for i, txt in enumerate(numeration):
    if txt!=6 and txt!=21 and txt!=23 :
        ax[1,0].annotate(txt, (logT[i]-0.02,logL[i]))
    ax[1,1].annotate(txt, (logT[i]-0.02,logL[i]))
    ax[1,2].annotate(txt, (logT[i]-0.02,logL[i]))
    hdu.writelines(star[i])
    hdu.writelines(' ')
    hdu.writelines(str(txt))
    hdu.writelines('\n')
hdu.close()
ax[1,0].plot([],[],' ',label=r'size: TESS, $\sigma_{Tp}$ [mmag]')
ax[1,1].plot([],[],' ',label=r'size: $H\alpha$, $\sigma_{\langle v^0\rangle}$')
ax[1,2].plot([],[],' ',label=r'size: Photospheric line, $\sigma_{\langle v\rangle}$ [km/s]')
ax[1,0].legend(framealpha=0.1,loc=4),ax[1,1].legend(framealpha=0.1,loc=4),ax[1,2].legend(framealpha=0.1,loc=4)
ax[1,0].set_xlim(4.73,3.95),ax[1,1].set_xlim(4.73,3.95),ax[1,2].set_xlim(4.73,3.95)



# In[mean_clear]
studied_stars='HD21291','HD34085','HD91316','HD199478'
V_alpha=[-200,190],[-200,200],[-300,300],[-400,300]
H_alpha=6562.80
plt.close('all')
plt.rcParams["figure.figsize"]=(8,6)
for star_lim in range(len(studied_stars)):
   star=studied_stars[star_lim]
   v_alpha_1,v_alpha_2=V_alpha[star_lim][0],V_alpha[star_lim][1]
   if os.path.exists('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/H_alpha_6562'):
       path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/H_alpha_6562'
   elif os.path.exists('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/H_alpha'):
       path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/H_alpha'
   files=os.listdir(path)
   common_fm=[]
   zm,HJD_zm=[],[]
   if os.path.exists(path+'/uniform_ZM.txt')==True:
       os.remove(path+'/uniform_ZM.txt')
   hdu=open(path+'/uniform_ZM.txt','a')
   figzm,axzm=plt.subplots()
   X_int,Y_int=[],[]
   if len(files)>55:
       criteria=50/len(files)
   else: 
       criteria=1
   M=[[],[]]
   H_alpha=6562.80
   for file in files: 
       if file.startswith('HD') and file.endswith('.txt'):
           wv,v,f=reading_line_file(path+'/'+file)
           x_l,v_l,y_l=wv[v>v_alpha_1],v[v>v_alpha_1],f[v>v_alpha_1]
           x_l,y_l,v_l=x_l[v_l<v_alpha_2],y_l[v_l<v_alpha_2],v_l[v_l<v_alpha_2]
           llh=np.random.uniform(0,1)
           if llh<criteria:
               axzm.plot(v,f,alpha=0.5)
           X_int.append(np.around(v,3)),Y_int.append(f)
           M[0].append(min(v)),M[1].append(max(v))
   X_int,Y_int=np.array(X_int),np.array(Y_int)
   X_m,Y_m=[],[]
   for j in range(len(X_int)):
       Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
       X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
   X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)     
   axzm.plot(X_mean,Y_mean,'k',markersize=1)
   
   plt.tick_params(axis='y', labelsize=20)
   plt.tick_params(axis='x', labelsize=20)  
   hdu.close()
   axzm.grid()
   axzm.set_title(star+r', $H\alpha$',fontsize=20)
   axzm.set_xlabel('Velocity [km/s]',fontsize=20),axzm.set_ylabel('Normalized flux',fontsize=20)
   axzm.axvline(x=v_alpha_1,c='g'),axzm.axvline(x=v_alpha_2,c='g')
   plt.savefig('../stars/'+star+'/SUMMARY/'+star+'_wind_big.png')




#%%

star,folder,line,v1,v2='HD31327','4567',4567.84,-90,90
plt.rcParams["figure.figsize"]=(8,6)
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
files=os.listdir(path)
common_fm=[]
#line=line
fm,HJD_fm=[],[]
if os.path.exists(path+'/uniform_fm.txt')==True:
    os.remove(path+'/uniform_fm.txt')
hdu=open(path+'/uniform_fm.txt','a')
X_int,Y_int=[],[]
if len(files)>55:
    criteria=50/len(files)
else: 
    criteria=1
M=[[],[]]
for file in files: 
    if file.startswith('HD') and file.endswith('.txt'):
        wv,v,f=reading_line_file(path+'/'+file)
        x_l,v_l,y_l=wv[v>v1],v[v>v1],f[v>v1]
        x_l,y_l,v_l=x_l[v_l<v2],y_l[v_l<v2],v_l[v_l<v2]
        llh=np.random.uniform(0,1)
        X_int.append(np.around(v,3)),Y_int.append(np.around(f,3))
        M[0].append(np.around(min(v),3)),M[1].append(np.around(max(v),3))
        FM=first_moment(x_l,y_l,line)/zero_moment(x_l,y_l,line)
        fm.append(FM)
        hjd=file.split('_')[0].replace(star,'')
        hjd=float(hjd.replace('p','.'))
        HJD_fm.append(hjd)
        if 'SONG' in file:
            ins='SONG'
        if 'MERCATOR' in file:
            ins='MERCATOR'
        if 'FIES' in file:
            ins='FIES'
        info=FM,hjd,ins
        info=str(info)
        hdu.writelines(info)
        hdu.writelines('\n')
X_int,Y_int=np.array(X_int),np.array(Y_int)
X_m,Y_m=[],[]
for j in range(len(X_int)):
    Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
    X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)    

if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
    os.mkdir('../stars/'+star+'/SUMMARY')

hdu=open(path+'/uniform_fm.txt','r')
FM_old,HJD_old_fm,INS_old_fm=[],[],[]
for j in hdu:
    dat=j.split()
    FM_old.append(float(dat[0][1:-1]))
    HJD_old_fm.append(float(dat[1][:-1]))
    INS_old_fm.append(dat[2][:-1])
FM_old=np.array(FM_old)-np.mean(FM_old)
HJD_old_fm=np.array(HJD_old_fm)-2.45e6
INS_old_fm=np.array(INS_old_fm)
FM_i,HJD_fm_i=[[],[],[]],[[],[],[]]
FM_r,HJD_fm_r=[[],[],[]],[[],[],[]]
hdu.close()
for fm in range(len(FM_old)):
    if 'MERCATOR' in INS_old_fm[fm]:
        ins=0
    elif 'FIES' in INS_old_fm[fm]:
        ins=1
    elif 'SONG' in INS_old_fm[fm]:
        ins=2
    if abs(FM_old[fm])<2.5*np.std(FM_old):
        FM_i[ins].append(FM_old[fm]),HJD_fm_i[ins].append(HJD_old_fm[fm])
    else:
        FM_r[ins].append(FM_old[fm]),HJD_fm_r[ins].append(HJD_old_fm[fm])
FM,HJD_fm=np.concatenate(FM_i),np.concatenate(HJD_fm_i)        
np.random.seed(5)
plt.plot(HJD_fm,FM,'.')
plt.grid()

fig,axs=plt.subplots()
HJD_fm_1,FM_1=HJD_fm[HJD_fm<8922],FM[HJD_fm<8922]
HJD_fm_1,FM_1=HJD_fm_1[HJD_fm_1>8766],FM_1[HJD_fm_1>8766]
axs.plot(HJD_fm_1,FM_1,'.')
axs.grid()
axs.set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$'),axs[0].set_xlabel('HJD-2450000')
print(np.std(FM_1))
#%%
star='HD30614'
line=None
l1,l2=5888.5,5891
x='wv'
title='Insterstellar lines'
spectrum=findstar(star)
step=0.015
files=[ins for ins in spectrum if '_N_' in ins]
# Interpolate a function in the same wavelengths all of them. 
try:
    del X_int,Y_int,HJD,RV,snr
except NameError:
    pass
c= 299792.458 
X_int,Y_int,HJD,RV,snr=[],[],[],[],[]
plt.figure()
bar = progressbar.ProgressBar(max_value=len(files))
I=0
step=0.015
for s in files:
    #rv0=RV0(line,s,verbose=False,line=line)
    
    try:
        SM=sp.spec(s,line=line)
        if l1==None: 
            l1=SM.wave[0]
        if l2==None: 
            l2=SM.wave[-1]
        SM.resamp(step,lwl=l1,rwl=l2)
        wv,fl=SM.wave,SM.flux
        if x=='wv' or x=='wave' or x=='wavelength':
            plt.plot(wv,fl)
        elif x=='vel' or x=='velocity':
            v_x=(np.array(wv)-line)*c/line 
            plt.plot(v_x,fl)
    except: 
        pass
    I=I+1
    bar.update(I)    
#plt.xlim(l1,l2), plt.ylim(0,2)  
plt.grid()
if x=='wv' or x=='wave' or x=='wavelength':
    plt.xlabel(r'Wavelength [$\AA$]',fontsize=16),plt.ylabel(r'Normalized flux',fontsize=16)
elif x=='vel' or x=='velocity':
    plt.xlabel(r'Velocity [km/s]'),plt.ylabel(r'Normalized flux')
if title==None:
    plt.title(star+', '+instrument+', '+str(line))
else:
    plt.title(str(title),fontsize=16)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(axis='x', labelsize=16)  
    
    
    
    
    
    