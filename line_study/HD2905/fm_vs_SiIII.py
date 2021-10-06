#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 17:11:11 2021

@author: charlie
"""
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
from line_study import *





def info(file):
    f1 = open(file,'r')  
    HJD1_m,fm1_m,error_fm1_m=[],[],[]
    HJD1_f,fm1_f,error_fm1_f=[],[],[]
    
    for lin in f1: 
        l = lin.split() 
        if l[2].replace(',','')=="'MERCATOR'":
            HJD1_m.append(float(l[1].replace(',',''))),fm1_m.append(float(l[9].replace(',','')))
            error_fm1_m.append(float(l[10].replace(',','')))
        elif l[2].replace(',','')=="'FIES'":
            HJD1_f.append(float(l[1].replace(',',''))),fm1_f.append(float(l[9].replace(',','')))
            error_fm1_f.append(float(l[10].replace(',','')))
    lis=sorted(zip(HJD1_m,fm1_m,error_fm1_m))
    HJD1_m,fm1_m,error_fm1_m=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]
    lis=sorted(zip(HJD1_f,fm1_f,error_fm1_f))
    HJD1_f,fm1_f,error_fm1_f=[i[0] for i in lis],[i[1] for i in lis],[i[2] for i in lis]
    
    return HJD1_m,HJD1_f,fm1_m,fm1_f,error_fm1_m,error_fm1_f
# In[]
SiIII_1='../stars/HD2905/4552/00_HD2905_4552.txt'  
SiIII_2='../stars/HD2905/5739/00_HD2905_5739.txt'
H_alpha='../stars/HD2905/H_alpha_6562/00_HD2905_H_alpha_6562.txt'
H_beta='../stars/HD2905/H_beta_4861/00_HD2905_H_beta_4861.txt'
HeI='../stars/HD2905/HeI_5875/00_HD2905_HeI_5875.txt'


HJD1_m,HJD1_f,fm1_m,fm1_f,error_fm1_m,error_fm1_f=info(SiIII_1)
HJD_m_SiIII,HJD_f_SiIII,fm_m_SiIII,fm_f_SiIII,error_fm_m_SiIII,error_fm_f_SiIII=info(SiIII_2)    
   
HJD_m_alpha,HJD_f_alpha,fm_m_alpha,fm_f_alpha,error_fm_m_alpha,error_fm_f_alpha=info(H_alpha)
HJD_m_beta,HJD_f_beta,fm_m_beta,fm_f_beta,error_fm_m_beta,error_fm_f_beta=info(H_beta)  
HJD_m_HeI,HJD_f_HeI,fm_m_HeI,fm_f_HeI,error_fm_m_HeI,error_fm_f_HeI=info(HeI)  

plt.close('all')

plt.figure()
plt.errorbar(fm1_m,fm_m_SiIII,xerr=error_fm1_m,yerr=error_fm_m_SiIII,fmt='.',label='MERCATOR Si III 5739')
plt.errorbar(fm1_f,fm_f_SiIII,xerr=error_fm1_f,yerr=error_fm_f_SiIII,fmt='.',label='FIES Si III 5739')
plt.errorbar(0.9727661740853849,0.1725828482534277,xerr=0.5538682407529805,yerr=1.3911488934052363,fmt='.',label='Extreme')
plt.grid()
plt.xlabel('Si III (4552)'),plt.ylabel('First moment')
plt.legend()

plt.figure()
plt.errorbar(fm1_m,fm_m_alpha,xerr=error_fm1_m,yerr=error_fm_m_alpha,fmt='.',label=r'MERCATOR H$\alpha$')
plt.errorbar(fm1_f,fm_f_alpha,xerr=error_fm1_f,yerr=error_fm_f_alpha,fmt='.',label=r'FIES H$\alpha$')
plt.grid()
plt.xlabel('Si III (4552)'),plt.ylabel('First moment')
plt.legend()

plt.figure()
plt.errorbar(fm1_m,fm_m_beta,xerr=error_fm1_m,yerr=error_fm_m_beta,fmt='.',label=r'MERCATOR H$\beta$')
plt.errorbar(fm1_f,fm_f_beta,xerr=error_fm1_f,yerr=error_fm_f_beta,fmt='.',label=r'FIES H$\beta$')
plt.grid()
plt.xlabel('Si III (4552)'),plt.ylabel('First moment')
plt.legend()

plt.figure()
plt.errorbar(fm1_m,fm_m_HeI,xerr=error_fm1_m,yerr=error_fm_m_HeI,fmt='.',label=r'MERCATOR He I')
plt.errorbar(fm1_f,fm_f_HeI,xerr=error_fm1_f,yerr=error_fm_f_HeI,fmt='.',label=r'FIES He I')
plt.grid()
plt.xlabel('Si III (4552)'),plt.ylabel('First moment')
plt.legend()


#plt.savefig('./stars/HD2905/4552_5739.png')

plt.figure()
plt.errorbar(fm_m_alpha,fm_m_beta,xerr=error_fm_m_alpha,yerr=error_fm_m_beta,fmt='.',label=r'MERCATOR')
plt.errorbar(fm_f_alpha,fm_f_beta,xerr=error_fm_f_alpha,yerr=error_fm_f_beta,fmt='.',label=r'FIES')
plt.grid()
plt.xlabel(r'H$\alpha$'),plt.ylabel(r'H$\beta$')
plt.legend()
# In[] 

# Estudio de un caso extremo de FIES para ver si tuviera algo raro en el perfil
# Si III 5739.

arch='../stars/HD2905/5739/HD29052457594p44634_FIES_5739p734.txt'

#HD29052457594p44634_FIES_5739p734.txt
def get_data(file):    
    fdatos = open(file,'r')
    wvl,vlc,flx = [],[],[]
    for lin in fdatos: 
        l = lin.split()
        wvl.append(float(l[0])),vlc.append(float(l[1])),flx.append(float(l[2]))
    return wvl,vlc,flx
wave,vel,flux=get_data(arch)
plt.figure()
plt.plot(wave,flux)
plt.xlabel(r'$\lambda(\AA)$'),plt.ylabel(r'Normalized flux')
# Vemos que lo que tiene es una línea... telúrica(?)
