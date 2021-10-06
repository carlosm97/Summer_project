# In[header]
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 10:49:06 2021

@author: charlie
"""
import astropy.io.fits as pf
import numpy as np 
import matplotlib.pyplot as plt
import os
os.chdir('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/pyIACOB')  
from scipy.optimize import curve_fit
from scipy import integrate
import formulas as fmr
from scipy.interpolate import interp1d
from db import *
from spec import *
from RV import *
from line_study import *
import time
import progressbar
# Let's play with SONG's files... 

files=os.listdir('../HD2905_SONG')

# In[Pruebas]

#El 3 o el 4 es la longitud de onda, aunque desconozco qué es el otro... a hablar 
#con Sergio. Por el ejemplo, usaremos el 3 (TRES) como wavelength solution before the science exposure 

#El 2 es el blaze.


# Plot order 30 (counting stats at 0) with the wavelength solution before the science exposure
#plt.close('all')
plt.figure()
plt.plot(data[0,13,:])
plt.grid()

plt.figure()
plt.plot(data[1,13,:])
plt.grid()

plt.figure()
plt.plot(data[2,13,:])
plt.grid()

plt.figure()
plt.plot(data[3,13,:],'.')

plt.plot(data[4,13,:],'.')
plt.grid()

#Las keywords para el HJD y la VHelio son
#BJD-MID =
#BVC     =
plt.close('all')
fig1,ax1=plt.subplots(nrows=1, ncols=1)
fig3,ax3=plt.subplots(nrows=1, ncols=1)

cc_wv,cc_flux=data[3,0,:],data[0,0,:]/data[2,0,:]
ax3.plot(data[3,0,:],data[1,0,:]/data[0,0,:],'r')
for i in range(1,np.shape(data)[1]):
    cc_wv,cc_flux=np.concatenate((cc_wv,data[3,i,:])),np.concatenate((cc_flux,data[0,i,:]/data[2,i,:]))
    #ax1.plot(data[3,i,:],data[0,i,:]/data[2,i,:],'b')
    #ax2.plot(data[3,i,:],data[1,i,:]/data[2,i,:],'r')
    ax3.plot(data[3,i,:],data[1,i,:]/data[0,i,:],'r')
ax1.plot(cc_wv,cc_flux)
ax1.grid(),ax2.grid()#,ax3.grid()


'''
# Este sería Si III 4552.622
plt.figure()
plt.plot(data[3,4,:],data[0,4,:]/data[2,4,:])
plt.plot(data[3,5,:],data[0,5,:]/data[2,5,:])

# Para ver a qué orden corresponde algo 
plt.figure()
for i in range(0,np.shape(data)[1]):
    plt.plot(data[3,i,:],np.linspace(i,i,np.shape(data)[2]))
'''

# In[]

# Vamos a tratar de normalizar un orden particular "a lo bruto", a falta de hablar 
# con Sergio para ver diferencias entre las 2 o 3 cosas vistas. 

# Por similitud con el ejemplo, miramos 1 y 3 (los mismos usados para concatenar)

# Orden 13, correspondiente a H alpha

plt.close('all')
fig,ax=plt.subplots()
wave,flux=data[3,13,:],data[0,13,:]/data[2,13,:]

flux_cont,wave_cont = np.concatenate((flux[wave<4855],flux[wave>4867])),np.concatenate((wave[wave<4855],wave[wave>4867]))
    
plt.figure()
plt.plot(wave,flux)
c_fit = np.poly1d(np.polyfit(wave_cont,flux_cont, 2))
plt.plot(wave,c_fit(wave))
plt.grid()
flux_norm=flux/c_fit(wave)
ax.plot(wave,flux_norm)
fn=np.concatenate((flux_norm[wave<4855],flux_norm[wave>4867]))

#plt.figure()
#plt.plot(wave_cont,fn)
sigma=np.std(fn)
fn=flux_norm*(flux_norm>1-sigma)*(flux_norm<1+sigma)
wave_copy=wave.copy()
wave_copy[fn==0.]=np.nan
fn[fn==0.]=np.nan
ax.plot(wave_copy,np.linspace(1.+3*sigma,1.+3*sigma,len(wave_copy)),'.',markersize=3)
ax.grid()


# In[H BETA]
'''H BETA'''
# Usamos 1 en lugar de 0 porque con 0 había algunos que no existían, y como íbamos 
# a normalizar igual...
#Vamos a intentar jugar a algo: 
c=299792.458 
'''OJO! Hay que sumar la heliocéntrica para corregir bien...'''
plt.close('all')
fig,ax=plt.subplots()
fig1,ax1=plt.subplots()
j=0
z=0
STD,fil=[],[]
WV,FLN,BJD=[],[],[]
bar = progressbar.ProgressBar(max_value=len(files))
for i in files:
    name = '../HD2905_SONG/'+i
    data = pf.getdata( name )
    hdr = pf.getheader(name)
    odr = order(data,hdr,4780)#6562.8)
    wave,flux=data[3,odr,:],data[0,odr,:]/data[2,odr,:]
    if (sum(abs(flux))==0)==True:
        flux=data[1,odr,:]/data[2,odr,:]
    if np.mean(flux)/np.std(flux)>10: 
        fn,cont,std=norm_SONG(wave,flux)
        STD.append(std),fil.append(i)
        ax1.plot(wave*(1+hdr.get('BVC')/c),fn)        
        odr=order(data,hdr,4862)
        wave,flux=data[3,odr,:],data[0,odr,:]/data[2,odr,:]
        if (sum(abs(flux))==0)==True:
            flux=data[1,odr,:]/data[2,odr,:]
        #print(i)
        fn,cont,std=norm_SONG(wave,flux)
        STD.append(std),fil.append(i)
        ax.plot(wave*(1+hdr.get('BVC')/c),fn)
        WV.append(wave),FLN.append(fn)
        BJD.append(hdr.get('BJD-MID'))
    else:
        z+=1
    j+=1
    bar.update(j)
    #print(j,'-',len(files))
print(z)
ax.grid()
ax.set_xlabel(r'$\lambda[\AA]$'),ax.set_ylabel('Normalized flux')
ax.set_xlim(4840,4880),ax.set_ylim(0.,2.)



# In[Caso ruidoso]
# Caso muy ruidoso: i=fil[794]
fig1,ax1=plt.subplots()

i=fil[794]
name = '../HD2905_SONG/'+i
data = pf.getdata( name )
if len(data[3,13,:][data[3,13,:]<4860])!=0:
    wave,flux=data[3,13,:],data[1,13,:]/data[2,13,:]
elif len(data[3,12,:][data[3,12,:]<4860])!=0:
    wave,flux=data[3,12,:],data[1,12,:]/data[2,12,:]
fn,cont,std=norm_SONG(wave,flux)
ax1.plot(wave,fn)
ax1.grid()
ax1.set_xlabel(r'$\lambda[\AA]$'),ax1.set_ylabel('Normalized flux')
ax1.set_xlim(4840,4880),ax1.set_ylim(0.6,1.8)

# In[Prueba SiIII]
# 4552.622 cae en el extremo//conexión de 2 órdenes.
plt.close('all')
i=files[0]
name = '../HD2905_SONG/'+i
data = pf.getdata( name )
hdr = pf.getheader(name)
#wave,flux=data[3,5,:],data[0,5,:]/data[2,5,:]


fn,cont,std,cnt=norm_SONG(wave,flux,continuum=True)
plt.plot(wave,fn)

odr=order(data,hdr,5739)


wave,flux=data[3,odr,:],data[0,odr,:]/data[2,odr,:]
fn,cont,std=norm_SONG(wave,flux)
plt.plot(wave,fn)
plt.grid()

plt.figure()

i=files[0]
name = '../HD2905_SONG/'+i
data = pf.getdata( name )

wave,flux=data[3,33,:],data[0,33,:]/data[2,33,:]


fn,cont,std,cnt=norm_SONG(wave,flux,continuum=True)
plt.plot(wave,fn)
plt.grid()

''' # Comprobación de que la heliocéntrica debe sumarse. 
plt.close('all')
fig1,ax1=plt.subplots()
j=0
STD,fil=[],[]
viejo,nuevo=[],[]
for i in files:
    name = '../HD2905_SONG/'+i
    data = pf.getdata(name)
    hdr = pf.getheader(name)
    if len(data[3,35,:][data[3,35,:]>5890])!=0:
        wave,flux=data[3,35,:],data[1,35,:]/data[2,35,:]
        nuevo.append(int(i[3:7]))
    else: #if len(data[3,34,:][data[3,34,:]<5890])!=0:
        wave,flux=data[3,36,:],data[1,36,:]/data[2,36,:]
        viejo.append(int(i[3:7]))
    #else: #len(data[3,15,:][data[3,15,:]<4860])==0:
        #wave,flux=data[3,15,:],data[1,15,:]/data[2,15,:]
    #    print(i)
    #print(i)
    fn,cont,std=norm_SONG(wave,flux)
    STD.append(std),fil.append(i)
    ax1.plot(wave*(1+hdr.get('BVC')/c),fn)
    #ax.plot(wave,cont,'.')
    j+=1
    print(j,'-',len(files))
ax1.grid()
ax1.set_xlabel(r'$\lambda[\AA]$'),ax1.set_ylabel('Normalized flux')
ax1.set_ylim(-0.5,0.1),ax2.set_xlim(5888,5898)
'''
# In[SiIII]
from line_study import *
# Let's study Si III.

plt.close('all')
#fig,ax=plt.subplots()
#fig1,ax1=plt.subplots()
j,z=0,0
STD,fil=[],[]
for i in files:
    name = '../HD2905_SONG/'+i
    data = pf.getdata( name )
    hdr = pf.getheader(name)
    ordr = order(data,hdr,4780)#6562.8)
    wave,flux=data[3,ordr,:],data[0,ordr,:]/data[2,ordr,:]
    if (sum(abs(flux))==0)==True:
        flux=data[1,ordr,:]/data[2,ordr,:]
    if np.mean(flux)/np.std(flux)>15:
        fn,cont,std=norm_SONG(wave,flux)
        #STD.append(np.std(flux)),fil.append(i)
        #ax1.plot(wave*(1+hdr.get('BVC')/c),fn)  
        ordr = order(data,hdr,4567.840)
        wave,flux=data[3,ordr,:],data[0,ordr,:]/data[2,ordr,:]
        if (sum(abs(flux))==0)==True:
            flux=data[1,ordr,:]/data[2,ordr,:]
        try: 
            fn,cont,std=norm_SONG(wave,flux,interactive=True)    
            STD.append(std)
        except TypeError:
            pass
        if 1/std>20:
            #ax.plot(wave*(1+hdr.get('BVC')/c),fn)
            pass
        else:
            j+=1
            print(i)
        #break 
    else:
        j+=1
    z+=1
    print(z,'-',len(files))
print(j)
#ax.grid()
#ax.set_xlabel(r'$\lambda[\AA]$'),ax.set_ylabel('Normalized flux')
#ax.set_ylim(0.,2.)


'''
# Este espectro resulta muy ruidoso, así que lo apartamos (por el momento)
i='s1_2016-08-10T04-19-09_ext.fits'
fig,ax=plt.subplots()
fig1,ax1=plt.subplots()
name = '../HD2905_SONG/'+i
data = pf.getdata( name )
hdr = pf.getheader(name)
ordr = order(data,hdr,4780)#6562.8)
wave,flux=data[3,ordr,:],data[1,ordr,:]/data[2,ordr,:]
if np.mean(flux)/np.std(flux)>21.5 and np.mean(flux)/np.std(flux)<21.6:
    fn,cont,std=norm_SONG(wave,flux)
    #STD.append(np.std(flux)),fil.append(i)
    ax1.plot(wave*(1+hdr.get('BVC')/c),fn)  
    ordr = order(data,hdr,4567.840)
    wave,flux=data[3,ordr,:],data[1,ordr,:]/data[2,ordr,:]
    fn,cont,std=norm_SONG(wave,flux)
    STD.append(std)
    if 1/std>20:
        ax.plot(wave*(1+hdr.get('BVC')/c),fn,label=i)
        ax.plot(wv_cnt,np.linspace(1.3,1.3,len(wv_cnt)),'.')
'''


# In[]

# Añadir en ese punto algo que, si la línea está muy al borde, muestre un aviso
# y pregunte si seleccionar esa o no. 

# Vamos a tratar de crear archivos en formato fits equivalentes a los FIES o MERCATOR. 







# In[continuacion]  
''' Continuamos con lo guardado. Yo, H beta.'''
c=299792.458 
HJD=[]
fn_common=[]
line=4861.325 
lim1,lim2=1.05,1.2
c1,c2=250,-200
wave_common=np.arange(4850,4870,0.015)
plt.close('all')
for i in range(len(WV)):
    inter=interp1d(WV[i],FLN[i]) 
    y=inter(wave_common)
    pol=np.polyfit(np.concatenate((wave_common[:c1],wave_common[c2:])),np.concatenate((y[:c1],y[c2:])),1)
    corr=y/np.polyval(pol,wave_common)
    if 1/np.nanstd(np.concatenate((corr[:c1],corr[c2:])))>10: # We avoid data with really bad SNR
        fn_common.append(corr)
        HJD.append(BJD[i])

X_int,Y_int=wave_common,np.vstack(fn_common)
X_mean,Y_mean=wave_common,np.mean(Y_int,axis=0)
Error=np.nanstd(Y_int,axis=0)
v_x=(X_mean-line)*c/line 
plt.figure()
plt.errorbar(v_x,Y_mean,yerr=Error,fmt='.',color='k',label='Mean spectrum')
[plt.plot(v_x,Y_int[i,:],alpha=0.3) for i in range(0,len(X_int))]
plt.xlabel('Velocity [km/s]'),plt.ylabel('Normalized flux'),plt.title(line)
plt.grid()
plt.ylim(0.7,1.1),plt.xlim(v_x[0],v_x[-1])
plt.legend()

# Computing the TSV and deciding sigma. 
std1=[np.nanstd(Y_int[i][:c1]) for i in range(0,len(Y_int))]
std2=[np.nanstd(Y_int[i][c2:]) for i in range(0,len(Y_int))]
STD=(np.array(std1)+np.array(std2))/2
std=np.mean(STD)

plt.figure()
plt.plot(v_x,Error)
plt.axhline(y=lim2*np.mean(std),color='g',label=str(lim2)+' mean std')
plt.axhline(y=lim1*np.mean(std),color='r',label='Cut criteria ('+str(lim1)+' mean std)')
plt.ylabel('TVS')
plt.xlabel('Velocity [km/s]')
plt.legend()
plt.grid()
xg,yg=[],[]
M = int(len(Error)/2)
for i in range(0,M):
    mmax_1p1=X_mean[M+i]
    if Error[M+i]>lim2*np.mean(std):
        continue
    else: 
        break
for i in range(0,M):
    mmin_1p1=X_mean[M-i]
    if Error[M-i]>lim2*np.mean(std):
        continue
    else: 
        break
    
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
# Vamos a escoger un criterio para sigma tal que la FWHM de la distribución sea 
# justamente la correspondiente a 1. 
sigma=((mmin_1p1-mmin_1)+(mmax_1-mmax_1p1))/2/2.355
# Pasarlo a velocidades para plotear

plt.axvline(x=(mmin_1-sigma*2.355-line)*c/line,color='g')
plt.axvline(x=(mmin_1+sigma*2.355-line)*c/line,color='g')

plt.axvline(x=(mmax_1-sigma*2.355-line)*c/line,color='g')
plt.axvline(x=(mmax_1+sigma*2.355-line)*c/line,color='g')
plt.title(str(line))

# The study of the line itself.

first_m,error_fm=[],[]
rv,error_rv=[],[]
I=0
BJD_list=[]
bar = progressbar.ProgressBar(max_value=len(Y_int))
for j in range(len(Y_int)):
    v,mini=[],[]
    #BJD,std,x_Si,y_Si=fmr.get_data('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_prueba/'+j)
    #BJD_list.append(BJD)
    #Z=sorted(zip(BJD_list,m1_def,mini))
    #BJD_list,m1_def,mini=[Z[i][0] for i in range(0,len(BJD_list))],[Z[i][1] for i in range(0,len(BJD_list))],[Z[i][2] for i in range(0,len(BJD_list))]
    x_Si,y_Si= X_int, Y_int[j]
    for i in range(0,50):
        rmmin=np.random.normal(loc=mmin_1,scale=2*sigma)
        rmmax=np.random.normal(loc=mmax_1,scale=2*sigma)
        ran1,ran2=fmr.find_nearest(x_Si,rmmin)[1],fmr.find_nearest(x_Si,rmmax)[1]
        m0=zero_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],line)
        m1=first_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],line)
        v.append(m1/m0)
        try:
            mini.append(radial_velocity(x_Si[ran1:ran2],y_Si[ran1:ran2],line,[0.2,5,0]))
        except RuntimeError:
            print('RuntimeError for ',x_Si[ran1],x_Si[ran2])
            continue
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
    bar.update(I)
plt.figure()
plt.errorbar(np.array(HJD),first_m,yerr=error_fm,fmt='.',label='First moment')
plt.errorbar(np.array(HJD),rv,yerr=error_rv,fmt='.',label='Gaussian radial velocity')
#plt.plot(np.array(HJD)-2450000,RV,'.',label='Program radial velocity')
plt.grid()
plt.tick_params(axis='x', labelsize=17),plt.tick_params(axis='y', labelsize=17)
plt.ylabel('Velocity',fontsize=17), plt.xlabel(r'BJD-2450000 [days]',fontsize=17)
plt.legend()
'''
plt.close('all')
plt.figure()'''
#desv_rv=np.std(np.array(rv)-np.array(RV)) # Generamos un error tal que ambas medidas
# sean coherentes
'''
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
#plt.close('all')'''
#if np.mean(desv_rv)>=0.5*np.mean(rv):
#    desv_rv=error_rv
plt.figure()
#plt.errorbar(first_m_g,rv_g,xerr=error_fm_g,yerr=error_rv_g,fmt='.',label='Gaussian')
mask=abs(np.array(first_m)-np.array(rv))<3*np.nanmean(abs(np.array(first_m)-np.array(rv)))
plt.errorbar(first_m,rv,xerr=error_fm,yerr=error_rv,fmt='.',label='TVS')
plt.xlabel('First moment [km/s]')
plt.ylabel('Radial velocity [km/s]')
plt.plot([min(rv),max(rv)],[min(rv),max(rv)],label='1-1 line')
plt.grid()   

linear_model=np.polyfit(np.array(first_m)[mask],np.array(rv)[mask],1)
linear_model_fn=np.poly1d(linear_model)
pendiente=linear_model[0]
D_pendiente=np.sqrt(len(first_m))*max(error_fm)/np.sqrt(len(first_m)*sum(np.array(first_m)**2)-(sum(first_m))**2)
print('Slope=',pendiente,r'$\pm$',D_pendiente)
plt.plot(first_m,linear_model_fn(first_m),color="green", label='Fit line')
plt.title('Slope='+str(round(pendiente,3))+r'$\pm$'+str(round(D_pendiente,3))+r', $\lambda$'+str(line))
plt.legend()


#%%
runcell('[header]', '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/SONG.py')
print('Header done')
runcell('[H BETA]', '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/SONG.py')
print('H done')
runcell('[continuacion]', '/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/SONG.py')
print('Line study done')
#%%
import tarfile
import os
import shutil
from db import *
# Fast code to extract tgz folders. 
path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/Spectra/'
list_dir=os.listdir(path)
'''
for ld in list_dir: 
    if ld.endswith('tgz'):
        try:
            file = tarfile.open(path+ld)
            file.extractall('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/Spectra')
            file.close()
            original=path+ld
            target=path+'OLD/'+ld
            shutil.move(original,target)
        except:
            print('Review ',ld)
'''
# Star list+avilable spectrum number for each. 
s_lc=open('../line_study/stars.txt','r')
stars_lc=[]
for s in s_lc:
    dat=s.split()
    stars_lc.append(s.replace('\n',''))
s_lc.close()
sv,sl=[],[]
try: os.remove('../line_study/star_list.txt')
except: pass
star_list=open('../line_study/star_list.txt','a')
for ld in list_dir: 
    if ld.startswith('HD'):
        star=ld[:ld.find('_')]
        if star not in sv and star in stars_lc:
            sv.append(star)
            sl.append((len(findstar(star)),star))
sl=sorted(sl,reverse=True)

for star in sl:
    star_list.writelines(star[1])
    star_list.writelines('   ')
    star_list.writelines(str(star[0]))
    star_list.writelines('\n')
star_list.close()
#star_list=sorted(star_list,reverse=True)

study_stars='HD36861','HD91316','HD188209','HD34078','HD2905'
for s in study_stars:
    print(np.vstack(sl)[np.where(np.vstack(sl)==s)[0]])












