#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 15:58:57 2021

@author: charlie
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
os.chdir('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/pyIACOB')
from scipy.optimize import curve_fit

Ej=1

def lectura(fichero):
    dat=[]
    fich_1 = open(fichero)
    for lin in fich_1:
        l = lin.split()
        dat.append(l)
    fich_1.close
    return(dat)   

def select(Dat,cantidad1):
    A = Dat[0]
    Object1= A.index(cantidad1)
    def Observable(objeto):
        res = []
        for i in range(3,403):
            res.append(float(Dat[i][objeto]))
        return(res)   
    Obser1 = Observable(Object1)
    return(Obser1)
path='/home/charlie/Desktop/MSc/Primero/Primer cuatrimetre/Estructuras y evolucion estelar/Eassy/Modelos_Sergio'
modelos = os.listdir(path+'/modelos/')
modelos.sort()
# If you don't remember the name of any object you want to plot, uncomment the
# next lines to check an example. 
plt.close('all')   
   
ZAMS=[[],[]]
fig, axs = plt.subplots()
axs.invert_xaxis()

M0,LOG_L=[],[]
logTc,logrhoc=[],[]
for mod in modelos:
    if mod[-5]=='0': #or mod[-5]=='4':
        if mod[1:4]=='005' or mod[1:4]=='009' or  mod[1:4]=='015' or mod[1:4]=='025'\
        or mod[1:4]=='040'or mod[1:4]=='060' or mod[1:4]=='085':
            Fich=lectura(path+'/modelos/'+mod)
            Log_L,Log_Teff=select(Fich,'lg(L)'),select(Fich,'lg(Teff)')
            ZAMS[0].append(Log_Teff[0])
            ZAMS[1].append(Log_L[0]-np.log10(float(str(mod[2:4]))))
            axs.plot(Log_Teff[:190],Log_L[:190]-np.log10(float(str(mod[2:4]))),c='k')
            axs.text(Log_Teff[110]+0.03,Log_L[110]-np.log10(float(str(mod[2:4]))),str(mod[2:4])+r'$M_{\odot}$')
            #plt.xlabel(r'log($T_{eff}$)')
            #plt.ylabel('log(L)')
            #plt.legend()
            #plt.title('H-R with initial rotation velocity 0')
ZAMS[0].sort(),ZAMS[1].sort()
axs.plot(ZAMS[0],ZAMS[1],'k')
axs.grid()
axs.set_xlabel(r'log($T_{eff}$[K])'),axs.set_ylabel(r'log($L/L_{\odot}$)')
axs.set_xlim(4.7,4),axs.set_ylim(3.1,4.3)
#axs.legend(fontsize=18)
hdu = open('../stars_Carlos.txt')
ID, Teff, lgf =[],[],[]
for j in hdu:
    i=j.split()
    if i[0].startswith('HD'):
        ID.append(i[0]),Teff.append(float(i[1])),lgf.append(float(i[2]))
        axs.plot(4+np.log10(float(i[1])),5.39-(float(i[2])),'.',label=i[0],markersize=10)
    
#axs.plot(4+np.log10(Teff),5.39-np.vstack(lgf),'.',label=ID)
axs.legend()
#plt.xlabel('log(Teff)'),plt.ylabel('L_spec')
hdu.close()



plt.show()













