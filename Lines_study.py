#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: charlie
"""
import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
import math as mth
import os
from scipy.optimize import curve_fit
from scipy import integrate
plt.close('all')
plt.rcParams['figure.figsize'] = [8, 6]

c = 299792.458
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx
def ajuste_mas(fit):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    y=hdul[0].data.copy()
    first = head['CRVAL1']
    step = head['CDELT1']
    vel = head['I-VBAR']
    beta = vel/c
    x2 = np.linspace(3763.75/(1+beta),(3763.75+len(y[0])*0.015625)/(1+beta),len(y[0]))
    hdul.close()
    return x2, y[0]
def ajuste_menos(fit):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    y=hdul[0].data.copy()
    first = head['CRVAL1']
    step = head['CDELT1']
    vel = head['I-VBAR']
    beta = vel/c
    x2 = np.linspace(3763.75/(1-beta),(3763.75+len(y[0])*0.015625)/(1-beta),len(y[0]))
    hdul.close()
    return x2, y[0]
def gaussian_flux(x, a, sigma, mu):
    # Function to create a flux as a gaussian. It it normalized at 1. 
    r = 1+a/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))
    return(r)  
def adjust(function, data, pars0):
    # Function to adjust 3 values (a, sigma, mu). To do a more clear code. 
    pars, cov = curve_fit(f=function, xdata=data[0], ydata=data[1], p0=pars0, bounds=(-np.inf, np.inf))
    return(pars, cov) 
def Pars0(dat):   
    # Let's calculate here a first approximation for the initial values to iniciate 
    # the fit. It could improve the computing time.
    Mu = sum(dat[0])/len(dat[0])
    sigma2 = sum(dat[0]**2) /len(dat[0]) - Mu**2
    Sig = sigma2**0.5   
    Max = max(dat[1]-1)*Sig*np.sqrt(2*np.pi)    
    pars0 = [Max, Sig, Mu]
    return(pars0)
def fit_gaussian(central_ojo,ristra):
    x,y=ristra[0],ristra[1]
    E=find_nearest(x,central_ojo)
    A = max(abs(1-y[E[1]-100:E[1]+100]))
    E_min=find_nearest(y,min((y[E[1]-100:E[1]+100])))
    E_max=find_nearest(y,max((y[E[1]-100:E[1]+100])))
    if (E_max[0]-1)>(1-E_min[0]):
        E_min=E_max
    mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))
    y_E=y[E_min[1]-400:E_min[1]+400]*mask
    x_E=x[E_min[1]-400:E_min[1]+400]*mask
    
    # Ajuste de gaussiana al perfil 
    X_E=[x_E[x_E !=0],y_E[y_E !=0]]
    pairs, covar=adjust(gaussian_flux,X_E,Pars0(X_E))
    x_lin=np.linspace(x_E[x_E !=0][0],x_E[x_E !=0][-1],100)
    return [x_lin,gaussian_flux(x_lin,*pairs),pairs]

def elem(fich,table):
    fichero= open(fich,'r')
    for lin in fichero:
        l = lin.split()
        try:
            table.append(float(l[0]))
        except:
            np.nan
    fichero.close()
    return table
def get_info(fit,key):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    info=head[key]
    del head
    hdul.close()
    return info
def minimo(central_ojo,ristra):
    x,y=ristra[0],ristra[1]
    E=find_nearest(x,central_ojo)
    A = max(abs(1-y[E[1]-100:E[1]+100]))
    E_min=find_nearest(y,min((y[E[1]-100:E[1]+100])))
    E_max=find_nearest(y,max((y[E[1]-100:E[1]+100])))
    if (E_max[0]-1)>(1-E_min[0]):
        E_min=E_max
    mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))#y[E_min[1]-400:E_min[1]+400]<=0.93
    return E_min

def linea(fich,x,y,l1_1,l1_2,l2_1,l2_2,central_ojo,destino):
    # l1_1,l1_2,l2_1,l2_2 son los límites del continuo.
    # Escogemos las zonas del continuo para renormalizarlas:  4546:4548 y 4556:4558 en el caso Si III
    c1=[find_nearest(x,l1_1)[1],find_nearest(x,l1_2)[1]]
    c2=[find_nearest(x,l2_1)[1],find_nearest(x,l2_2)[1]]
    x_c1,x_c2,y_c1,y_c2=x[c1[0]:c1[1]],x[c2[0]:c2[1]],y[c1[0]:c1[1]],y[c2[0]:c2[1]]
    x_c,y_c=np.concatenate((x_c1,x_c2)),np.concatenate((y_c1,y_c2))
    P=np.polyfit(x_c,y_c,2)
    Y_fit=np.polyval(P, x)
    x,y,Y_fit=x[c1[0]:c2[1]],y[c1[0]:c2[1]],Y_fit[c1[0]:c2[1]]
    info=[]
    std=np.std(y_c)
    snr=1/std
    BJD=get_info(fich,'BJD')
    # Para escoger la línea, cogemos un criterio de la media de 10 ptos. con media alejada de 1 más de 1 std.
    X=[x,y/Y_fit]
    MIN=minimo(central_ojo,X)
    y_gr,x_gr,y_lo,x_lo=[],[],[],[]
    for j in range(0,len(X[0])):
        x_gr.append(X[0][MIN[1]+j]),y_gr.append(X[1][MIN[1]+j])
        if np.mean(X[1][MIN[1]+j-5:MIN[1]+j+5])>1-std:
            break
    for j in range(1,len(X[0])):
        x_lo.append(X[0][MIN[1]-j]),y_lo.append(X[1][MIN[1]-j])
        if np.mean(X[1][MIN[1]-j-5:MIN[1]-j+5])>1-std:
            break
    x_line,y_line=np.concatenate((np.array((x_lo[::-1])),np.array(x_gr))),\
    np.concatenate((np.array((y_lo[::-1])),np.array(y_gr)))
    # Esta es la información que queremos guardar de cada línea
    info.append(x_line),info.append(y_line),info.append(std),info.append(BJD)
    
    BJD_s=str(BJD)
    BJD_s=BJD_s.replace('.','p')
    np.savetxt(destino+BJD_s+'MERCATOR_SiIII'+str(4552)+'.txt',(info),fmt='%s') 
def zero_moment(x_l,y_l,line):
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int0=integrate.cumtrapz(F_l, v_l, initial=0)    
    return y_int0[-1]    

def first_moment(x_l,y_l,line):
    BJD,std,x_l,y_l=get_data('./SiIII/'+i)
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int1 = integrate.cumtrapz(F_l*v_l, v_l, initial=0)
    return y_int1[-1]

def radial_velocity(x_l,y_l,line):
    v_l = (x_l-line)*c/line    
    return fit_gaussian(0,[v_l,y_l])[-1][-1]

def second_moment(x_l,y_l,line,centre=None):
    BJD,std,x_l,y_l=get_data('./SiIII/'+i)
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line-centre
    F_l=1-y_l
    y_int2 = integrate.cumtrapz(F_l*v_l**2, v_l, initial=0)
    return y_int2[-1]
def third_moment(x_l,y_l,line,centre=None):
    BJD,std,x_l,y_l=get_data('./SiIII/'+i)
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line-centre
    F_l=1-y_l
    y_int3 = integrate.cumtrapz(F_l*v_l**3, v_l, initial=0)
    return y_int3[-1]
# In[]
# Creamos los ficheros con los datos que tenemos. Al hacerlo por función, es mucho más rápido.
plt.close('all')
contenido = os.listdir('./HD2905_FH')
M_2905 = []

for file in contenido:
    if file.endswith('M_V85000.fits'):
        M_2905.append(file)
for i in M_2905:
    if i=='HD2905_20110116_205853_M_V85000.fits':
        res = ajuste_mas('./HD2905_FH/'+i)
    else:
        res = ajuste_menos('./HD2905_FH/'+i)
    x1,x2=find_nearest(res[0],4545),find_nearest(res[0],4560)
    
    x,y=res[0][x1[1]:x2[1]],res[1][x1[1]:x2[1]]
    linea('./HD2905_FH/'+i,x,y,4546,4548,4556,4558,4552,'./SiIII_prueba/HD_2905')    
# In[]
    
''' 
plt.close('all')
contenido = os.listdir('./SiIII_prueba')
M_2905_SiIII = []
# Sacamos los datos de la línea de todos los ficheros. 
for file in contenido:
    M_2905_SiIII.append(file)
f_menos=[]
BJD_list = []
for i in M_2905_SiIII:
    BJD,std,x_l,y_l=get_data('./SiIII_prueba/'+i)
    globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    plt.plot(globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII'])
    BJD_list.append(BJD)
plt.grid()
plt.tick_params(axis='x', labelsize=17)
plt.tick_params(axis='y', labelsize=17)
plt.xlabel(r'Wavelength [$\AA$]',fontsize=17),
plt.ylabel('Normalized flux',fontsize=17)'''
# In[]
# Integramos y ploteamos.     
    
contenido = os.listdir('./SiIII')
M_2905_SiIII = []
# Sacamos los datos de la línea de todos los ficheros. 
for file in contenido:
    M_2905_SiIII.append(file)
    
    
BJD_list,m1_def,mini,m2_def,m3_def = [],[],[],[],[]
for i in M_2905_SiIII:
    SiIII=4552.622
    BJD,std,x_Si,y_Si=get_data('./SiIII/'+i)
    f0=zero_moment(x_Si,y_Si,SiIII)
    f1_norm=first_moment(x_Si,y_Si,SiIII)/f0
    f2_norm=second_moment(x_Si,y_Si,SiIII,f1_norm)/f0
    f3_norm=third_moment(x_Si,y_Si,SiIII,f1_norm)/f0
    BJD_list.append(BJD)
    m1_def.append(f1_norm),m2_def.append(f2_norm),m3_def.append(f3_norm)
    mini.append(radial_velocity(x_Si,y_Si,SiIII))
# Ordenamos las listas por fecha. 
Z=sorted(zip(BJD_list,m1_def,mini))
BJD_list,m1_def,mini=[Z[i][0] for i in range(0,len(BJD_list))],[Z[i][1] for i in range(0,len(BJD_list))],[Z[i][2] for i in range(0,len(BJD_list))]

    
    
plt.plot(np.array(BJD_list)-2450000,m1_def,'.',label='First moment',markersize=8)
plt.plot(np.array(BJD_list)-2450000,mini,'x',label='Radial velocity')
plt.grid()
plt.ylabel(r'Velocity [km/s]',fontsize=17),plt.xlabel(r'BJD-2450000 [days] ',fontsize=17)
plt.tick_params(axis='x', labelsize=17)
plt.tick_params(axis='y', labelsize=17)
plt.legend(fontsize=17)
plt.savefig('./image/SiIII_fm.png')

plt.figure()
plt.title('Intrinsic velocity (First moment-Radial velocity)',fontsize=17)
plt.plot(np.array(BJD_list)-2450000,np.array(m1_def)-np.array(mini),'.',markersize=8)
plt.grid(),plt.ylim(-15,15)
plt.ylabel(r'Velocity [km/s]',fontsize=17),plt.xlabel(r'BJD-2450000 [days] ',fontsize=17)
plt.tick_params(axis='x', labelsize=17)
plt.tick_params(axis='y', labelsize=17)
plt.savefig('./image/SiIII_fm_rv.png')


plt.figure()
plt.plot(np.array(BJD_list)-2450000,np.array(m1_def)-np.array(mini),'.',markersize=10)
plt.grid()
plt.ylabel(r'Velocity [km/s]',fontsize=17),plt.xlabel(r'BJD-2450000 [days] ',fontsize=17)
plt.tick_params(axis='x', labelsize=17)
plt.tick_params(axis='y', labelsize=17)
plt.title('Intrinsic velocity (First moment-Radial velocity)',fontsize=17)
plt.xlim(6640,6650),plt.ylim(-15,10)
plt.savefig('./image/SiIII_zoom_fm_rv.png')
np.mean(np.array(m1_def)-np.array(mini))

plt.figure()
plt.plot(np.array(BJD_list)-2450000,np.array(m1_def),'.',label='First moment',markersize=10)
plt.grid()
plt.ylabel(r'Velocity [km/s]',fontsize=17),plt.xlabel(r'BJD-2450000 [days] ',fontsize=17)
plt.tick_params(axis='x', labelsize=17)
plt.tick_params(axis='y', labelsize=17)
plt.plot(np.array(BJD_list)-2450000,np.array(mini),'x',label='Radial Velocity',markersize=10)
#plt.title('Intrinsic velocity (First moment-Radial velocity)',fontsize=17)
plt.xlim(6640,6650),plt.ylim(-20,10)
plt.legend(fontsize=17)
plt.savefig('./image/SiIII_zoom_fm.png')
#plt.axhline(y=np.mean(np.array(mini)),color='g')    
    
    
    
    
    