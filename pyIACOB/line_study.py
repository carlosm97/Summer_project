#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:03:10 2021

@author: charlie
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:26:36 2021

@author: charlie
"""
import os
os.chdir('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/pyIACOB')  
import time
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import integrate
import formulas as fmr
from scipy.interpolate import interp1d
from db import *
#from spec import *
import spec as sp
from RV import *
import progressbar
import time
c=299792.458 
def colored(r, g, b, text):
    '''Print in terminal in other color'''
    return "\033[38;2;{};{};{}m{} \033[38;2;255;255;255m".format(r, g, b, text)



def fit_gaussian(central_ojo,ristra,pars0=None):
    '''Fit a gaussian. 
    
    Parameters
    ----------
    central_ojo : float,
        Point close to the central value to start iteratting the fit. 
    
    ristra : list or array,
        Array of arrays containing x (ristra[0]) and y (ristra[1]) arrays to be fitted. 
     
    pars0 : array, optional
        Initial values to compute [Max, Sig, Mu]. If not, compute them with function Pars0
        
    Returns
    -------
    x_lin : array, 
        Value of x fitted
    
    y_lin : array,
        Value of y fitted. 
        
    pairs : array,
        Values of the peak, sigma and mu (in this order) from the fit. 
    '''
    x,y=ristra[0],ristra[1]
    #mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))#y[E_min[1]-400:E_min[1]+400]<=0.93    
    # Ajuste de gaussiana al perfil 
    X_E=[x,y]#[x_E[x_E !=0],y_E[y_E !=0]]
    if pars0==None:
        pars0=fmr.Pars0(X_E)
    
    pairs, covar=fmr.adjust(fmr.gaussian_flux,X_E,pars0)
    x_lin=x#np.linspace(x_E[x_E !=0][0],x_E[x_E !=0][-1],100)
    return [x_lin,fmr.gaussian_flux(x_lin,*pairs),pairs]
def radial_velocity(x_l,y_l,line,pars0=None):
    '''Apply fit_gaussian to get the radial velocity 
    
    Parameters
    ----------
    x_l : array,
        Wavelength array
        
    y_l : array,
        Normalized flux array
        
    line : float
        Wavelength in the rest frame of reference.
        
    pars0 : array, optional
        Initial values to compute [Max, Sig, Mu]. If not, compute them with function Pars0'
    
    Returns
    -------
    See 'fit_gaussian' help. 
        '''
    v_l = (x_l-line)*c/line    
    pars0[2]=(pars0[2]-line)*c/line  
    return fit_gaussian(0,[v_l,y_l],pars0)[-1][-1]  
def zero_moment(x_l,y_l,line):
    '''
    Parameters
    ----------
    x_l : array,
        Wavelength.

    y_l : array,
        Normalized flux.

    line : float,
        Centre of the line to be studied.

    Returns
    -------
    int : float,
        Numerical value of the integral
    '''
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int0=integrate.cumtrapz(F_l, v_l, initial=0)    
    return y_int0[-1]    

def zero_reading(star,folder):
    '''
    Function to read the files containing the zero moment where it is computed. 
    
    Parameters
    ----------
    
    star : str, 
        Name of the star to study
        
    folder : str,
        Name of the folder to of the line which zero moment want to be read. 
        
    Return
    ------
    
    OF : list, 
        List of the names of the original files of the spectrum
        
    LF : list,
        List of the names of the line files used to compute the momentum 
    
    HJD : array,
        Array containing the Julian data of the observations
        
    ZM : array,
        Array of the zero moment 
    '''
    path='../stars/'+star+'/'+str(folder)+'/00_zero_moment.txt'
    hdu=open(path)
    OF,LF,HJD,ZM=[],[],[],[]
    for l in hdu: 
        dat=l.split()
        OF.append(dat[0]),LF.append(dat[1]),HJD.append(float(dat[2].replace(',',''))),ZM.append(float(dat[3].replace(')','')))
    hdu.close()
    return OF,LF,np.vstack(HJD),np.vstack(ZM)



def first_moment(x_l,y_l,line):
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int1 = integrate.cumtrapz(F_l*v_l, v_l, initial=0)
    return y_int1[-1]
def order(dat,hdr,line):
    orde=None
    for i in range(np.shape(dat)[1]):
        if dat[3,i,0]*(1+hdr.get('BVC')/c)<line and  dat[3,i,-1]*(1+hdr.get('BVC')/c)>line :
            try:
                if min(line-dat[3,i,0]*(1+hdr.get('BVC')/c),dat[3,i,-1]*(1+hdr.get('BVC')/c)-line)>min(abs(line-dat[3,orde,0]*(1+hdr.get('BVC')/c)),abs(dat[3,orde,-1]*(1+hdr.get('BVC')/c)-line)):
                    orde=i
            except ValueError:
                orde=i
    return orde
def norm_SONG(wave,flux,continuum=False,interactive=False,criteria=0.05):
    cont=True
    flux_norm=flux
    pol=np.array([1,1,1])
    i=0
    sigma,sigma1=1,0
    while abs(sigma-sigma1)>1e-3*sigma:
        flux_cont,wave_cont=flux_norm*cont,wave*cont
        flux_cont[flux_cont==0]=np.nan
        wave_cont[wave_cont==0]=np.nan
        if len(flux_cont[np.isfinite(flux_cont)])<criteria*len(flux):
            if interactive==True:
                plt.close('Confirmation')
                plt.figure('Confirmation')
                plt.plot(wave,flux_norm)
                plt.plot(wave_cont,flux_cont,'.')
                plt.grid()
                plt.show()
                plt.pause(1)
                continuity=input(colored(255,0,0,'WARNING: Too many rejected points ('+str(len(flux)-len(flux_cont[np.isfinite(flux_cont)]))+'/'+str(len(flux))+')')+
                                 "\n" +'Do you want to continue normalizing? [discard/y/N] ')
                if continuity=='yes' or continuity=='Yes'  or continuity=='YES'  or continuity=='y'  or continuity=='Y':
                    continue
                elif continuity=='discard' or continuity=='Discard' or continuity=='DISCARD' or continuity=='d' or continuity=='D':
                    return 
                else:
                    break
            else:
                pass
        std=np.nanstd(flux_cont)
        idx = np.isfinite(flux_cont) & np.isfinite(wave_cont)
        pol=np.polyfit(wave_cont[idx],flux_cont[idx], 2)
        c_fit = np.poly1d(pol)
        try:    
            flux_norm=flux_norm/c_fit(wave)
        except TypeError:
            flux_norm=flux/c_fit(wave)
        sigma1=sigma
        sigma=np.std(flux_norm)
        #plt.plot(wave,flux_norm,label=i)
        cont=(flux_norm>1-std) & (flux_norm<1+std)
    flux_cont,wave_cont=flux_norm*cont,wave*cont
    flux_cont[flux_cont==0]=np.nan
    wave_cont[wave_cont==0]=np.nan
    std=np.nanstd(flux_cont)
    idx = np.isfinite(flux_cont) & np.isfinite(wave_cont)
    pol=np.polyfit(wave_cont[idx],flux_cont[idx], 2)
    c_fit = np.poly1d(pol)
    flux_norm=flux_norm/c_fit(wave)
    y=np.linspace(1.+2*std,1.+2*std,len(wave))
    y[(wave*cont)==0]=np.nan  
    if continuum==True:
        cont=wave_cont
        return flux_norm,y,std,cont
    else: 
        return flux_norm,y,std
   
def mean_line(star,folder,line,instrument):
    '''
    Function to read line files and compute the mean line. 
    
    Parameters
    
    ----------
    star : str,
        Name of the studied star.
        
    folder : str,
        Name of the folder of the line to get the mean.
        
    line : float,
        Line wavelenght at the rest frame of reference.
        
    instrument : str, list or tuple,
        Name of the instrument (or instruments) whose data have to be averaged.

    Returns
    -------
    WV : array,
        Array of wavelegths.
        
    V : array,
        Array of velocities (wavelength transformed).
        
    F : array,
        Array of used flux.
        
    HJD : array,
        Julian date of the observation.
        
    flux_mean : array,
        Array of mean flux.

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
        elif instrument==None and l.endswith('txt'):
            str_line=str(line).replace('.','p')
            l=l.replace('_FIES'+'_'+str_line+'.txt','')
            l=l.replace('_SONG'+'_'+str_line+'.txt','')
            l=l.replace('_MERCATOR'+'_'+str_line+'.txt','')
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
def prueba_linea(star,instrument,line,l1,l2,x='wv',title=None):
    '''
    Give a first idea of the parameters to be used in line study function
    
    Parameters
    ----------
    star : str,
        Name of the star to be studied.
        
    instrument : str,
        Name of the instrument which data are going to be studied
    
    line : float
        Central line to be studied 
        
    l1 & l2 : float
        Limits os the line to be studied 
    
    x: string, optional
        Unit in which we want to see the x axis. It can be seen in wavelength
        (def) or in velocity (also accepted 'vel').
    
    title: string, optional
        Title to be plotted in the figure. If none is given, the title is the 
        star+line+instrument. 
    
    Return
    ------
    Nothing. Generates figures of the line. 
    '''
    spectrum=findstar(star)
    try:
        del instrum
    except NameError:
        pass  
    if instrument=='MERCATOR' or instrument=='mercator' or instrument=='hermes' or instrument=='HERMES':
        condition='_M_'
        ext='_FH'
    elif instrument=='FIES' or instrument=='fies':
        condition='_N_'  
        ext='_FH'      
    elif instrument=='SONG' or instrument=='song':
        condition='s1_'  
        ext='_SONG'
    else: 
        print('Unknown instrument')
    files=[ins for ins in spectrum if condition in ins]
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
        plt.xlabel(r'Wavelength [$\AA$]'),plt.ylabel(r'Normalized flux')
    elif x=='vel' or x=='velocity':
        plt.xlabel(r'Velocity [km/s]'),plt.ylabel(r'Normalized flux')
    if title==None:
        plt.title(star+', '+instrument+', '+str(line))
    else:
        plt.title(str(title))
    return    
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
def line_study(star,instrument,line,ll,rl,step,lim1=1,lim2=1.2,MC=500,plt_mean=False,\
               target=None,c1=60,c2=-60,snr_min=0,save=True,automatic_selection=True,\
               ll1=None,ll2=None,sel='vel',telluric=False,sig=None,normalization=True,\
               spectrum=None):
    """
    Function to study the line and create files with the information.
    
    Parameters
    ----------
    star : str, 
        Name of the star to be studied. The same as used in function findstar.
    
    instrument : str,
        Name of the instrument the spectra are taken with 
    
    line : float,
        Lab wavelength of the line to be studied (in the used units)
    
    ll : float,
        Subjective left limit of the line
    
    rl : float,
        Subjective right limit of the line
    
    step: float,
        Step to interpolate all the spectra
    
    lim1 & lim2 : floats, optional
        Limits of the TVS selection criteria. Default is 1 and 1.2 std.
    
    MC : int, optional
        Number of iterations in the MC simulation. Default is 500.
    
    plt_mean : bool, optional 
        Plots all spectra line and the mean line in the same graph if True.
        Default is False
    
    target : str, optional
        Name of the folder to save the data obtained. If not, 
        it uses a folder whose name is the 4 first digits of the line. 
        
    c1 & c2: int, optional
        Size considered continuum from the beggining to c1 and from c2 to the end. 
    
    snr_min: float, optional
        Minimum signal to noise ratio to study and save the line in a given spectra.
        Default is 0.
        
    save: bool, optional
        Variable to choose if save the results or not. Default is true. Usefull 
        to run several times when looking for the best parameters.
        
    telluric: bool, optional 
        To correct telluric lines. Default is False. 
        
    AUTOMATIC SELECTION Y ESAS COSAS
        
    sig: float, optional 
        Manually select sigma value to stimate the uncertainty. It is selected 
        in the units of sel (default in wavelength). It can be changed to 'vel'.
        
    normalization: bool, optional
        If True (def), it normalizes the extremes of the line using the intervals
        :c1 at left and c2: at right. 
        
    Returns
    -------
        
    First moment and its error by MC stimation
    
    Radial velocity and its error by MC+error by getting different velocities
    
    Wavelength and flux interpolated.
    """
    plt.rcParams["figure.figsize"] = (8,6)
    if target==None:
        target=str(line)[:4]
    if os.path.exists('../stars/'+star+'/'+target)==False:
        os.mkdir('../stars/'+star+'/'+target)
    else: 
        pass
    if spectrum==None:
        spectrum=findstar(star)
    try:
        del instrum
    except NameError:
        pass
    if instrument=='MERCATOR' or instrument=='mercator' or instrument=='hermes' or instrument=='HERMES':
        condition='_M_'
        ext='_FH'
    elif instrument=='FIES' or instrument=='fies':
        condition='_N_'  
        ext='_FH'  
    elif instrument=='SONG' or instrument=='song':
        condition='s1_'  
        ext='_SONG'
    else:
        print('\n','Unknown instrument.')
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
    ll,rl=line-step*int((line-ll)/step),line+step*int((rl-line)/step)
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
    
    # Computing the TSV and deciding sigma. 
    std1=[np.std(Y_int[i][:c1]) for i in range(0,len(X_int))]
    std2=[np.std(Y_int[i][c2:]) for i in range(0,len(X_int))]
    STD=(np.array(std1)+np.array(std2))/2
    std=(np.mean(std1)+np.mean(std2))/2
    
    plt.figure()
    plt.plot(v_x,Error)
    plt.axhline(y=lim2*np.mean(std),color='g',label=str(lim2)+' mean std')
    plt.axhline(y=lim1*np.mean(std),color='r',label='Cut criteria ('+str(lim1)+' mean std)')
    plt.ylabel('TVS')
    plt.xlabel('Velocity [km/s]')
    plt.legend()
    plt.grid()
    xg,yg=[],[]
    if automatic_selection==True:
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
        # Vamos a escoger un criterio para sigma tal que la FWHM de la distribución sea 
        # justamente la correspondiente a 1. 
        if sig==None:
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
            sigma=((mmin_1p1-mmin_1)+(mmax_1-mmax_1p1))/2/2.355
            # Pasarlo a velocidades para plotear
        else:
            if sel=='vel' or sel=='velocity':
                sigma=line*(sig/c)
            else: 
                sigma=sig
        plt.axvline(x=(mmin_1-sigma*2.355-line)*c/line,color='g')
        plt.axvline(x=(mmin_1+sigma*2.355-line)*c/line,color='g')
        
        plt.axvline(x=(mmax_1-sigma*2.355-line)*c/line,color='g')
        plt.axvline(x=(mmax_1+sigma*2.355-line)*c/line,color='g')
        plt.title(str(line))
        plt.savefig('../stars/'+star+'/'+target+'/TVS_criteria_'+instrument+'.png')
    else:
        if type(ll1)==float or type(ll1)==int and type(ll2)==float or type(ll2)==int:
            if sel=='vel' or sel=='velocity' or sel=='kms':
                mmin_1=line*ll1/c+line
                mmax_1=line*ll2/c+line
            elif sel=='wv' or sel=='wave' or sel=='wavelength' or sel=='lam':
                mmin_1=ll1
                mmax_1=ll2
            else: 
                print('Not valid selection')
                return
        else:
            print('Not valid selection limit')
            return
        if sig==None:
            sigma=((mmax_1-mmin_1))/40
            # Pasarlo a velocidades para plotear
        else:
            if sel=='vel' or sel=='velocity':
                sigma=line*(sig/c)
            else: 
                sigma=sig
    first_m,error_fm=[],[]
    rv,error_rv=[],[]
    I=0
    BJD_list=[]
    bar = progressbar.ProgressBar(max_value=len(Y_int))
    if plt_mean==True:
        plt.figure()
        plt.errorbar(v_x,Y_mean,yerr=Error,fmt='.',color='k',label='Mean spectrum')
        [plt.plot(v_x,Y_int[i,:],alpha=0.3) for i in range(0,len(X_int))]
        plt.xlabel('Velocity [km/s]'),plt.ylabel('Normalized flux'),plt.title(line)
        plt.axvline(x=(mmin_1-line)*c/line,color='g')
        plt.axvline(x=(mmax_1-line)*c/line,color='g')       
        plt.grid()
        plt.xlim(v_x[0],v_x[-1])
        plt.title(str(line))
        plt.axvline(x=(mmin_1-sigma*2.355-line)*c/line,color='g',alpha=0.05)
        plt.axvline(x=(mmin_1+sigma*2.355-line)*c/line,color='g',alpha=0.05)
        
        plt.axvline(x=(mmax_1-sigma*2.355-line)*c/line,color='g',alpha=0.05)
        plt.axvline(x=(mmax_1+sigma*2.355-line)*c/line,color='g',alpha=0.05)
        plt.savefig('../stars/'+star+'/'+target+'/mean_'+instrument+'.pdf')    
        
    
    print('\n','Error stimation and first moment computation')
    for j in range(0,len(X_int)):
        v,mini=[],[]
        #BJD,std,x_Si,y_Si=fmr.get_data('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/line_prueba/'+j)
        #BJD_list.append(BJD)
        #Z=sorted(zip(BJD_list,m1_def,mini))
        #BJD_list,m1_def,mini=[Z[i][0] for i in range(0,len(BJD_list))],[Z[i][1] for i in range(0,len(BJD_list))],[Z[i][2] for i in range(0,len(BJD_list))]
        x_Si,y_Si= X_int[j], Y_int[j]
        for i in range(0,MC):
            rmmin=np.random.normal(loc=mmin_1,scale=2*sigma)
            rmmax=np.random.normal(loc=mmax_1,scale=2*sigma)
            ran1,ran2=fmr.find_nearest(x_Si,rmmin)[1],fmr.find_nearest(x_Si,rmmax)[1]
            m0=zero_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],line)
            m1=first_moment(x_Si[ran1:ran2],y_Si[ran1:ran2],line)
            v.append(m1/m0)
            '''
            try:
                mini.append(radial_velocity(x_Si[ran1:ran2],y_Si[ran1:ran2],line,[0.2,5,0]))
            except RuntimeError:
                print('RuntimeError for ',x_Si[ran1],x_Si[ran2])
                continue
            '''
        #plt.figure()     
        n, bins = np.histogram(v, bins=100,density='true')
        mids = 0.5*(bins[1:] + bins[:-1])
        mean = np.average(mids, weights=n)
        var = np.average((mids - mean)**2, weights=n)
        error = np.sqrt(var)
        first_m.append(mean),error_fm.append(error)   
        '''
        n, bins = np.histogram(mini, bins=100,density='true')
        mids = 0.5*(bins[1:] + bins[:-1])
        mean = np.average(mids, weights=n)
        var = np.average((mids - mean)**2, weights=n)
        error = np.sqrt(var)
        rv.append(mean),error_rv.append(error)
        '''
        ran1,ran2=fmr.find_nearest(x_Si,mmin_1)[1],fmr.find_nearest(x_Si,mmax_1)[1]
        try: rv.append(radial_velocity(x_Si[ran1:ran2],y_Si[ran1:ran2],line,[max(abs(1-y_Si[ran1:ran2])),5,float(np.array(x_Si[ran1:ran2])[np.array(abs(1-y_Si[ran1:ran2]))==max(abs(1-y_Si[ran1:ran2]))])]))
        except: rv.append(np.nan)
        #rv.append(radial_velocity(x_Si[ran1:ran2],y_Si[ran1:ran2],line,[0.2,5,0]))
        error_rv.append(np.nan)
        I=I+1
        bar.update(I)
    
    #plt.figure()
    #plt.errorbar(np.array(HJD)-2450000,first_m,yerr=error_fm,fmt='.',label='First moment')
    #plt.errorbar(np.array(HJD)-2450000,rv,yerr=error_rv,fmt='.',label='Gaussian radial velocity')
    #plt.plot(np.array(HJD)-2450000,RV,'.',label='Program radial velocity')
    #plt.grid()
    #plt.tick_params(axis='x', labelsize=17),plt.tick_params(axis='y', labelsize=17)
    #plt.ylabel('Velocity',fontsize=17), plt.xlabel(r'BJD-2450000 [days]',fontsize=17)
    #plt.legend()
    #
    #plt.close('all')
    #plt.figure()
    desv_rv=np.std(np.array(rv)-np.array(RV))
    
    if np.mean(desv_rv)>=0.5*np.mean(rv):
        desv_rv=error_rv
    plt.figure()
    #plt.errorbar(first_m_g,rv_g,xerr=error_fm_g,yerr=error_rv_g,fmt='.',label='Gaussian')
    plt.errorbar(first_m,rv,xerr=error_fm,yerr=desv_rv,fmt='.',label='TVS')
    plt.xlabel('First moment [km/s]')
    plt.ylabel('Radial velocity [km/s]')
    plt.plot([min(rv),max(rv)],[min(rv),max(rv)],label='1-1 line')
    plt.grid()   
    try:
        linear_model=np.polyfit(first_m,rv,1)
        linear_model_fn=np.poly1d(linear_model)
        pendiente=linear_model[0]
        D_pendiente=np.sqrt(len(first_m))*max(error_fm)/np.sqrt(len(first_m)*sum(np.array(first_m)**2)-(sum(first_m))**2)
        print('Slope=',pendiente,r'$\pm$',D_pendiente)
        plt.plot(first_m,linear_model_fn(first_m),color="green", label='Fit line')
        plt.title(str(line)+', Slope='+str(round(pendiente,3))+r'$\pm$'+str(round(D_pendiente,3))+r', $\lambda$'+str(line))
    except:
        plt.title(str(line))
        print('Linear fit do not converge')
    plt.legend()
    plt.savefig('../stars/'+star+'/'+target+'/fm_rv_'+instrument+'.png')
    plt.close('all')
    if save==True:
        # Guardamos las líneas en ficheros: 
        import csv
            
        for i in range(0,len(X_int)):
            info=zip(X_int[i],(X_int[i]-line)*c/line,Y_int[i])
            BJD_s=str(HJD[i])
            BJD_s=BJD_s.replace('.','p')
            line_s=str(line)
            line_s=line_s.replace('.','p')
            with open('../stars/'+star+'/'+target+'/'+star+BJD_s+'_'+instrument+'_'+line_s+'.txt', 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(info)
            f.close()
        # nombre, JD, intrumento, resolución, mínimo y máximo de la línea, sigma a su alrededor, 
        # 
        f=open('../stars/'+star+'/'+target+'/'+'00_'+star+'_'+target+'.txt', 'a')
        f.writelines('Original file,    Line file,     HJD,    Instrument,     Resolution,     lam_min,    lam_max,\
            step,   sigma,     snr,    first moment,   error fm,   Radial velocity,    error rv, \
            Other error rv')
        f.writelines('\n') 
        if instrument=='SONG' or instrument=='song':
            resolution=77000
        else:
            resolution=files[i][-10:-4]
        for i in range(0,len(X_int)):
            try:
                BJD_s=str(HJD[i])
                BJD_s=BJD_s.replace('.','p')
                info=files[i].replace('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/'\
                                  +star+ext+'/',''),\
                star+BJD_s+'_'+instrument+'_'+line_s+'.txt',\
                HJD[i],instrument,resolution,mmin_1,mmax_1,step,sigma,1/STD[i],first_m[i],error_fm[i],\
                rv[i],error_rv[i],desv_rv[i]
            except:
                info=files[i].replace('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/'\
                                  +star+ext+'/',''),\
                star+BJD_s+'_'+instrument+'_'+line_s+'.txt',\
                HJD[i],instrument,resolution,mmin_1,mmax_1,step,sigma,1/STD[i],first_m[i],error_fm[i],\
                rv[i],error_rv[i],desv_rv
            f.writelines(str(info))
            f.writelines('\n')
        f.close()
    return(first_m,error_fm,rv,error_rv,desv_rv,X_int,Y_int)    

class reading_summary():
    '''
    Read the summary of lines as a class.
    
    Parameters
    ----------
    star : str,
        Name of the star or path to the summary 
        
    folder : str, optional
        Name of the folder of the line to be read. Default is None. If star is 
        not a path to the line file and folder is not an input, the class ask 
        for it. 
        
    Return
    ------
        Information of the summary in lists of str or floats.
    '''
    def __init__(self, star,folder=None,verbose=True):
        try:
            hdu = open(star)
        except:
            if verbose==True:
                if folder==None:
                    folder = input('Please, write name of the folder: ')
                path = '../stars/'+star+'/'+str(folder)+'/00_'+star+'_'+str(folder)+'.txt'
                hdu = open(path)
            else: pass
        Original_file,line_file,HJD,Instrument,Resolution,lam_min,lam_max,step,sigma,snr,first_moment,\
        error_fm,Radial_velocity,error_rv,Other_error_rv=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        for j in hdu:
            try:
                dat=j.split()
                if len(dat)==15:
                    Original_file.append(dat[0].replace(',',''))
                    line_file.append(dat[1].replace(',',''))
                    HJD.append(float(dat[2].replace(',','')))
                    Instrument.append(dat[3].replace(',',''))
                    Resolution.append(float(dat[4].replace(',','').replace("'","")))
                    lam_min.append(float(dat[5].replace(',','')))
                    lam_max.append(float(dat[6].replace(',','')))
                    step.append(float(dat[7].replace(',','')))
                    sigma.append(float(dat[8].replace(',','')))
                    snr.append(float(dat[9].replace(',','')))
                    first_moment.append(float(dat[10].replace(',','')))
                    error_fm.append(float(dat[11].replace(',','')))
                    Radial_velocity.append(float(dat[12].replace(',','')))
                    error_rv.append(float(dat[13].replace(',','')))
                    Other_error_rv.append(float(dat[14].replace(',','')))
                if len(dat)==16:
                    Original_file.append(dat[0].replace(',',''))
                    line_file.append(dat[2].replace(',',''))
                    HJD.append(float(dat[3].replace(',','')))
                    Instrument.append(dat[4].replace(',',''))
                    Resolution.append(float(dat[5].replace(',','').replace("'","")))
                    lam_min.append(float(dat[6].replace(',','')))
                    lam_max.append(float(dat[7].replace(',','')))
                    step.append(float(dat[8].replace(',','')))
                    sigma.append(float(dat[9].replace(',','')))
                    snr.append(float(dat[10].replace(',','')))
                    first_moment.append(loat(dat[11].replace(',','')))
                    error_fm.append(float(dat[12].replace(',','')))
                    Radial_velocity.append(float(dat[13].replace(',','')))
                    error_rv.append(float(dat[14].replace(',','')))
                    Other_error_rv.append(float(dat[15].replace(',','')))
            except:
                continue
        self.Original_file=Original_file
        self.line_file=line_file
        self.HJD=HJD
        self.Instrument=Instrument
        self.Resolution=Resolution
        self.lam_min=lam_min
        self.lam_max=lam_max
        self.step=step
        self.sigma=sigma
        self.snr=snr
        self.first_moment=first_moment
        self.error_fm=error_fm
        self.Radial_velocity=Radial_velocity
        self.error_rv=error_rv
        self.Other_error_rv=Other_error_rv
def reading_line_file(file,path='',verbose=False):
    '''
    Parameters
    ----------
    file : TYPE
        DESCRIPTION.
    path : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    try:
        fichero= open(path+file,'r')
    except:
        if verbose==True:
            path=input('Please, write the path of the folder where the file is ')    
            fichero= open(path+file,'r')
        else:
            return None
    wv,v,f=[],[],[]
    for lin in fichero:
        l = lin.split()
        wv.append(float(l[0])),v.append(float(l[1])),f.append(float(l[2]))
    fichero.close()
    return np.array(wv),np.array(v),np.array(f)

def fm_fm(line2,line1='4552',star='HD2905',error=False):
    '''
    Function to get the comparison of different lines for HERMES, FIES, and SONG

    Parameters
    ----------
    
    line2 : str,
        Name of the folder of the second line to compare with.
        
    line1 : str, optional
        Name of the folder of the second line to compare with. The default is '4552'.
        
    star : str, optional
        Name of the star which lines wants to be compared. The default is 'HD2905'.
        
    error : bool, optional
        If true, it returns the error associated with the first moment. The default is False.

    Returns
    -------
    
    IF  error==False:
    correspondence_MERCATOR,correspondence_FIES,correspondence_SONG: matrix,
        3 matrix with the values of the correspondece for the different instrument. 
        If there are not data, matrix is empty. The first line is [:,0], while 
        the second is [:,1]
        
    IF error==True:
    correspondence_MERCATOR,error_MERCATOR,correspondence_FIES,error_FIES,correspondence_SONG,error_SONG    
        6 matrix. Added 3 of error with the same structure as the others. 
        
    '''
    correspondence_FIES,correspondence_MERCATOR,correspondence_SONG=[],[],[]
    error_FIES,error_MERCATOR,error_SONG=[],[],[]
    lin1=reading_summary(star,line1)
    lin2=reading_summary(star,line2)
    for fil in lin1.Original_file:
        try:
            i=np.where(np.array(lin2.Original_file)==fil)[0][0]
            o=np.where(np.array(lin1.Original_file)==fil)[0][0]
            if  lin1.Instrument[o]=="'FIES'":
                correspondence_FIES.append([lin1.first_moment[o],lin2.first_moment[i]])
                error_FIES.append([lin1.error_fm[o],lin2.error_fm[i]])
            elif lin1.Instrument[o]=="'MERCATOR'":
                correspondence_MERCATOR.append([lin1.first_moment[o],lin2.first_moment[i]])
                error_MERCATOR.append([lin1.error_fm[o],lin2.error_fm[i]])
            elif lin1.Instrument[o]=="'SONG'":
                correspondence_SONG.append([lin1.first_moment[o],lin2.first_moment[i]])
                error_SONG.append([lin1.error_fm[o],lin2.error_fm[i]])
        except:
            pass
    try:
        correspondence_SONG=np.vstack(correspondence_SONG)
        error_SONG=np.vstack(error_SONG)
    except:
        pass
    try:
        correspondence_MERCATOR=np.vstack(correspondence_MERCATOR)
        error_MERCATOR=np.vstack(error_MERCATOR)
    except:
        pass        
    try:
        correspondence_FIES=np.vstack(correspondence_FIES)
        error_FIES=np.vstack(error_FIES)        
    except:
        pass
    if error==False:
        return correspondence_MERCATOR,correspondence_FIES,correspondence_SONG       
    elif error==True:
        return correspondence_MERCATOR,error_MERCATOR,\
            correspondence_FIES,error_FIES,\
                correspondence_SONG,error_SONG       
  
def lectura(fichero):
    '''
    Function to read a file data. 
    
    Parameters
    ----------
    
    fichero: str,
        Name of the file to be read. 
        
    Returns
    -------
    
    dat: list, 
        List with the content of the file in str. 
    '''
    dat=[]
    fich_1 = open(fichero)
    for lin in fich_1:
        l = lin.split()
        dat.append(l)
    fich_1.close
    return(dat)   

def select(Dat,cantidad1):
    '''
    Function to select the data to be read of GINEBRA models

    Parameters
    ----------
    Dat : list,
        Result from function lectura
        
    cantidad1 : str
        quantity to be obtained

    Returns
    -------
    The information of the column selected.

    '''
    A = Dat[0]
    Object1= A.index(cantidad1)
    def Observable(objeto):
        res = []
        for i in range(3,403):
            res.append(float(Dat[i][objeto]))
        return(res)   
    Obser1 = Observable(Object1)
    return(Obser1)

path='../updated_northern_lc/new_lc_ascii'
lc=os.listdir(path)
def find_star_photometry(star):
    hdu = open('../line_study/star_list.txt','r')
    lst = [[],[],[]]
    for j in hdu:
        dat=j.split()
        lst[0].append(dat[0]),lst[1].append(dat[1]),lst[2].append(dat[2])
    hdu.close()
    return np.array(lst[2])[np.array(lst[0])==star][0]
def photometry(star):
    '''
    Function to get the photometric avilable data
    
    Parameters
    ----------
    
    star: str
        Name of the star to be found.
        
    Return 
    ------
    
    HJD: list,
        List of floats containing the observation data HJD-2.45e6
        
    mag: list,
        List of floats containing the measured magnitude (unit in mmag)
    '''
    
    for j in lc:
        if find_star_photometry(star) in j:
            st=j    
    path='../updated_northern_lc/new_lc_ascii'
    hdu = open(path+'/'+st)
    HJD,mag=[],[]
    for j in hdu:
        dat=j.split()
        if 'SC' in st or 'LC' in st or 'PDC' in st or 'SAP' in st: 
            HJD.append(float(dat[0])+7000.0),mag.append(float(dat[1]))
        elif 'ktwo' in st:
            HJD.append(float(dat[0])+4833.0),mag.append(float(dat[1]))
    hdu.close()
    return HJD,mag 
def gif(star,line_folder,date,F,v,mean,jd1,jd2,line_title=r'H$\alpha$, ',typ='Wind',ll1=None,ll2=None,pbt=None):
    '''
    Function to do a GIF with 6 plots: 3 rows and 2 columns.

    Parameters
    ----------
    
    star : str
        Name of the star to GIF. 
        
    line_folder : str
        Name of the folder of the line to do the GIF. 
        
    date : list
        List with the dates to plot.
        
    F : list
        List containing the flux for different dates. It could be the result of
        mean_line function.
        
    v : list
        List of velocities associated to the flux. It could be the result of
        mean_line function.
        
    mean : list
        List containing the mean line values. 
        
    jd1 & jd2: float
        Dates to zoom in. 
        
    line_title : TYPE, optional
        DESCRIPTION. The default is r'H$\alpha$'.
        
    typ : str, optional
        If 'wind' (default), the zero moment is used. If 'photosphere' or 
        'photospheric', it is used the first moment.

    ll1 & ll2 : float, optional
        Limits for the plotted TESS data. If None (default), it use the same 
        interval as for the line data. 

    Returns
    -------
    
    Nothing, but save the GIF. From top to bottom and left to right, the GIF 
    plots: the variation of the line for different dates, the photometric magnitude
    from TESS, the residues of the line for different dates divided by the mean line,
    zoom of the zero (defoult) or first moment, the HR diagram with the star and
    GINEBRA models, and the whole measured of the zero (defoult) or first moment.
    '''   
    
    plt.rcParams["figure.figsize"] = (8,6)
    Z=sorted(zip(date,F,v))
    DATE,FLUX=[Z[i][0] for i in range(0,len(date))],[Z[i][1] for i in range(0,len(date))]
    VELOCITY=[Z[i][2] for i in range(0,len(date))]
    SM=reading_summary(star,line_folder)
    ZM,JD=[],[]
    if os.path.exists('../stars/'+star+'/'+line_folder+'/uniform_ZM.txt')==True:
        opn=open('../stars/'+star+'/'+line_folder+'/uniform_ZM.txt','r')
        Z_HJD,Z_ZM=[],[]
        for l in opn:
            l=l.replace('(','')
            l=l.replace(')','')
            l=l.replace(',','')
            j=l.split()
            Z_HJD.append(float(j[1])),Z_ZM.append(float(j[0]))
        opn.close()
        print('\n exists \n')
    else:
        Z_OF,Z_LF,Z_HJD,Z_ZM=zero_reading(star,line_folder)
    Z_ZM=np.array(Z_ZM)
    Z_ZM=(Z_ZM-np.mean(Z_ZM))/np.mean(Z_ZM)
    for i in range(len(np.array(SM.HJD))):
        if Z_HJD[i]>jd1 and Z_HJD[i]<jd2:
            JD.append(Z_HJD[i]),ZM.append(Z_ZM[i])
    date,flux,velocity=[],[],[]
    Z=sorted(zip(JD,ZM))
    JD,ZM=[Z[i][0] for i in range(0,len(JD))],[Z[i][1] for i in range(0,len(JD))]
    for i in JD:
        if round(float(i),4) in np.around(DATE,4):
            if pbt!=None: 
                if pbt<np.random.rand():
                    continue
            pos=np.where(np.around(DATE,4)==round(float(i),4))[0][0]
            date.append(float(i)),flux.append(FLUX[pos])
            velocity.append(VELOCITY[pos])
    date=np.vstack(date)
    flux=np.vstack(flux)
    fecha=np.zeros(np.shape(velocity))
    for i in range(len(fecha[0,:])):
        for j in range(len(fecha[:,0])): 
            fecha[j,i]=date[j]
    import sys
    from matplotlib.animation import FuncAnimation
    try:
        plt.close(fig)
    except:
        pass
    fig,axes= plt.subplots(3,2)#,sharex=True)
    fig.subplots_adjust(right=0.85)
    #fig.set_tight_layout(True)
    ax1, ax5=axes[0]
    ax2, ax4=axes[1]
    ax6, ax3=axes[2]
    ax2.set_xlabel(r'Velocity [km/s]')
    ax1.set_ylabel(r'Normalized flux')
    ax2.set_ylabel('HJD [days]-2450000')
    cm = plt.cm.get_cmap('RdYlBu')
    im = ax2.scatter(velocity,fecha-2.45e6,c=flux/mean,s=0.5,cmap=cm,vmax=1.05)
    # Query the figure's on-screen size and DPI. Note that when saving the figure to
    # a file, we need to provide a DPI for that separately.
    print('fig size: {0} DPI, size in inches {1}'.format(
        fig.get_dpi(), fig.get_size_inches()))
    x = velocity[0] 
    line, = ax1.plot(x, flux[0], 'r-', linewidth=2)
    line2 = ax2.axhline(y=sorted(fecha[:,0])[0]-2.45e6,c='r')
    #ax1.sharex(ax2)
    ax1.grid()
    #ax1.set_ylim(0.8,1.6)
    ax1.plot(x,mean,'grey')
    ax1.plot(x,np.mean(flux,0),'k')

    #ax1.plot(x,np.mean((flux),0),'g')
    ax2.set_xlim(min(velocity[0]),max(velocity[0]))
    ax1.set_xlim(min(velocity[0]),max(velocity[0]))
    # Para compartir eje
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    if typ=='Wind' or typ=='wind':
        ax4.plot(np.array(JD)-2.45e6,ZM,'.')
        ax4.set_ylabel(r'$\Delta\langle v^0\rangle/\langle v^0\rangle_0$',labelpad=10)
        ax3.plot(np.array(Z_HJD)-2.45e6,Z_ZM,'.')
        ax3.set_ylabel(r'$\Delta\langle v^0\rangle/\langle v^0\rangle_0$',labelpad=0)
    elif typ=='Photosphere' or typ=='Photospheric' or typ=='photosphere' or typ=='photospheric':
        # First moment 
        lin=reading_summary(star,line_folder)
        FM=[]
        for f in date:
            FM.append(np.array(lin.first_moment)[np.array(lin.HJD)==f[0]])
        ax4.plot(np.array(JD)-2.45e6,FM,'.')
        ax4.set_ylabel(r'$\langle v\rangle$ [km/s]',labelpad=10)  
        ax3.plot(np.array(SM.HJD)-2.45e6,SM.first_moment,'.')
        ax3.set_ylabel(r'$\langle v\rangle$ [km/s]',labelpad=0)
    line3 = ax4.axvline(x=date[0][0]-2.45e6,c='r')
    line4 = ax5.axvline(x=date[0][0]-2.45e6,c='r')
   
    
    ax3.grid()
    ax3.set_xlabel('HJD [days]-2450000')
    
    ax3.axvline(x=sorted(fecha[:,0])[0]-2.45e6,c='g')
    ax3.axvline(x=sorted(fecha[:,0])[-1]-2.45e6,c='g')
    ax4.grid()
    
    #cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.95)
    fig.tight_layout(w_pad=6)
    cb_ax = fig.add_axes([0.45, 0.415, 0.02, 0.22]) #[left, bottom, width, height] 
    cbar = fig.colorbar(im, cax=cb_ax)
    #cbar.set_ticks(np.arange(0, 2, 0.5))
    #cbar.set_ticklabels(['low', 'medium', 'high'])
    HJD,mag=photometry(star)
    ax5.plot(HJD,np.array(mag)*1000,'.',markersize=0.5)
    ax5.invert_yaxis()
    ax5.grid(),ax5.set_ylabel(r'$\Delta$Tp [mmag]')
    if ll1==None and ll2==None:
        ax5.get_shared_x_axes().join(ax5, ax4)
        ll1,ll2=np.array(JD)[0]-2.45e6,np.array(JD)[-1]-2.45e6
        ax5.set_xticklabels([]) 
    ax5.set_xlim(ll1,ll2)
    ax4.set_xlim(jd1-2.45e6,jd2-2.45e6)
    #ax1.set_xscale('log')  
    #ax1.set_yscale('log')
    #ax1.set_ylim(1e-15,1e-1)
    path='/home/charlie/Desktop/MSc/Primero/Primer cuatrimetre/Estructuras y evolucion estelar/Eassy/Modelos_Sergio'
    modelos = os.listdir(path+'/modelos/')
    modelos.sort()
    # If you don't remember the name of any object you want to plot, uncomment the
    # next lines to check an example. 
    ZAMS=[[],[]]
    ax6.invert_xaxis()
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
                ax6.plot(Log_Teff[:200],log_L_g[:200],c='k')
                if float(mod[1:4])<80 and float(mod[1:4])>9:
                    ax6.text(4.2,log_L_g[190]+0.05,str(mod[2:4])+r'$M_{\odot}$')        
    ZAMS[0].sort(),ZAMS[1].sort()
    ax6.plot(ZAMS[0],ZAMS[1],'k')
    ax6.grid()
    ax6.set_xlabel(r'log($T_{eff}$[K])'),ax6.set_ylabel(r'log($L/L_{\odot}$)')
    ax6.set_xlim(4.73,4),ax6.set_ylim(3.1,4.3)
    #axs.legend(fontsize=18)
    hdu = open('../line_study/stars.txt')
    ID, Teff, lgf =[],[],[]
    for j in hdu:
        i=j.split()
        if i[0].startswith(star):
            ID.append(i[0]),Teff.append(float(i[1])),lgf.append(float(i[2]))
            ax6.plot(4+np.log10(float(i[1])),5.39-(float(i[2])),'.',label=i[0],markersize=10)
    #axs.plot(4+np.log10(Teff),5.39-np.vstack(lgf),'.',label=ID)
    #plt.xlabel('log(Teff)'),plt.ylabel('L_spec')
    hdu.close()
    bar = progressbar.ProgressBar(max_value=np.shape(flux)[0])
    def update(i):
        #label = 'shell {0}'.format(i+Kgen_8000[0][0])
        title = line_title+str(date[i][0].round(3))
        #print(label)
        # Update the line and the ax1es (with a new xlabel). Return a tuple of
        # "artists" that have to be redrawn for this frame.
        line.set_ydata(flux[i])
        line2.set_ydata(sorted(fecha[:,0])[i]-2.45e6)
        line3.set_xdata(date[i][0]-2.45e6)
        line4.set_xdata(date[i][0]-2.45e6)
        fig.suptitle(title)
        fig.subplots_adjust(top=0.95)
        bar.update(i)
        return line, ax1
    if __name__ == 'line_study' or __name__ == '__main__' :
        # FuncAnimation will call the 'update' function for each frame; here
        # animating over 10 frames, with an interval of 200ms between frames.
        anim = FuncAnimation(fig, update, frames=np.shape(flux)[0], interval=200)
        if len(sys.argv) > 1 and sys.argv[1] == 'save':
            anim.save('line.gif', dpi=80, writer='imagemagick')
        else:
            # plt.show() will just loop the animation forever.
            plt.show()
        anim.save('/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/GIFS/'+line_folder+'_'+star+'.gif', writer='imagemagick', fps=3)
        print('all')
    return None
  
def statistical_study(star,folder,line,v1,v2,v_alpha_1,v_alpha_2,ph_name='Si III',TESS=True):
   '''
   Function to do the statistical study of the star. It studies H_alpha for the 
   wind, and the specified line for the photosphere, as well as the photometry 
   
   Parameters
   ----------
   
   star : str, 
       Name of the star to be studied.
       
   folder : str, 
       Name of the folder of the line to be studied.
       
   line : float, 
       Wavelength of the line to be studied.
       
   v1 & v2 : floats,
       Commun limits (left and right) in velocity to integrate the photometric
       line and get a first moment without instrument differenciation.
   
   v_alpha_1,v_alpha_2 : floats,
       Commun limits (left and right) in velocity to integrate the wind line 
       and get a zero moment without instrument differenciation.
   
   ph_name : str, optional
       Element name of the studied line. Deffault is Si III
   TESS : bool, optional
       If TESS is True (def), it is assumed that there are TESS data, and are 
       ploted together with the others.
   Returns:
   -------
   
   sigma_fm : float,
       Sigma of the photospheric line.
   
   sigma_zm : float,
       Sigma of the wind line (H_alpha).
   mean_ZM: float
   
   sigma_ph : float,
       Sigma of the photometric data.
   len_fm
   len_zm    
   In addition to the previous information, it saves several figures: first moment
   vs HJD, and the study with different ammount of data
   '''
   
   
   plt.rcParams["figure.figsize"]=(8,6)
   path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star+'/'+folder
   files=os.listdir(path)
   common_fm=[]
   #line=line
   fm,HJD_fm=[],[]
   if os.path.exists(path+'/uniform_fm.txt')==True:
       os.remove(path+'/uniform_fm.txt')
   hdu=open(path+'/uniform_fm.txt','a')
   figfm,axfm=plt.subplots()
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
   axfm.plot(X_mean,Y_mean,'k',markersize=1)
   axfm.grid()
   axfm.set_title(star+r', $\lambda_0=$'+str(line)+r'$\AA$')
   axfm.set_xlabel('Velocity [km/s]',fontsize=16),axfm.set_ylabel('Normalized flux',fontsize=16)
   axfm.axvline(x=v1,c='g'),axfm.axvline(x=v2,c='g')
   axfm.tick_params(axis='y', labelsize=16),axfm.tick_params(axis='x', labelsize=16)  
   if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
       os.mkdir('../stars/'+star+'/SUMMARY')
   plt.savefig('../stars/'+star+'/SUMMARY/'+star+'_photosphere.png')

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
   #pp_fm_10,sigma_fm_10=[],[]
   pp_fm_25,sigma_fm_25=[],[]
   pp_fm_50,sigma_fm_50=[],[]
   for i in range(25000):
       # Selection of the data to be used:
       num=[]
       rnd=np.random.randint(0,len(HJD_fm),int(0.25*len(HJD_fm)))
       for j in rnd:
           num.append(float(FM[j]))
       #pp_fm_25.append(max(num)-min(num))
       sigma_fm_25.append(np.std(num))
       num=[]
       rnd=np.random.randint(0,len(HJD_fm),int(0.5*len(HJD_fm)))
       for j in rnd:
           num.append(float(FM[j]))
       #pp_fm_50.append(max(num)-min(num))
       sigma_fm_50.append(np.std(num))
   sigma_fm=np.std(FM)
   
   STD_1=[]
   step1,step10,step25,step50 = 1,2,7,50
   data_0_fm=min(HJD_fm)
   data_1_fm=max(HJD_fm)
   n_step_1,n_step_10,n_step_25=int((data_1_fm-data_0_fm)/step1)+1,int((data_1_fm-data_0_fm)/step10)+1,int((data_1_fm-data_0_fm)/step25)+1
   n_step_50=int((data_1_fm-data_0_fm)/step50)+1
   for i in range(n_step_1):
       mask1=np.logical_and(np.vstack(HJD_fm)>(i+data_0_fm),np.vstack(HJD_fm)<(i+step1+data_0_fm))
       if len(np.vstack(FM)[mask1])>=5:
           STD_1.append(np.std(np.vstack(FM)[mask1]))
   plt.close('all')
   sp,sa=plt.subplots(3)
   #print(len(STD_1),len(sigma_fm_25),len(sigma_fm_50))
   if len(STD_1)!=0:
       maxi,mini=max(max(STD_1),max(sigma_fm_25),max(sigma_fm_50)),min(min(STD_1),min(sigma_fm_25),min(sigma_fm_50))
   else:
       maxi,mini=max(max(sigma_fm_25),max(sigma_fm_50)),min(min(sigma_fm_25),min(sigma_fm_50))
   H_1, bins_1 = np.histogram(STD_1,bins=50,range=[mini,maxi])
   sa[0].bar((bins_1[1:]+bins_1[:-1])/2,H_1/max(H_1),alpha=0.3,width=np.mean(bins_1[1:]-bins_1[:-1]),label='1 day')
   sa[0].legend()
   H_f_25, bins_f_25 = np.histogram(sigma_fm_25,bins=50,range=[mini,maxi])
   H_m_25, bins_m_25 = np.histogram(sigma_fm_50,bins=50,range=[mini,maxi])
   sa[1].bar((bins_f_25[1:]+bins_f_25[:-1])/2,H_f_25/max(H_f_25),alpha=0.3,width=np.mean(bins_f_25[1:]-bins_f_25[:-1]),label='25%')
   sa[1].bar((bins_m_25[1:]+bins_m_25[:-1])/2,H_m_25/max(H_m_25),alpha=0.3,width=np.mean(bins_m_25[1:]-bins_m_25[:-1]),label='50%')
   sa[1].legend()
   sa[1].set_xlabel(r'$\sigma_{\langle v\rangle}$ [km/s]')
   sa[1].axvline(x=sigma_fm)
   sa[0].get_shared_x_axes().join(sa[0], sa[1])
   sa[2].plot(HJD_old_fm,FM_old,'.r',markersize=1)
   sa[2].plot(HJD_fm,FM,'.b',markersize=1)
   sa[2].grid()
   sa[2].set_xlabel(r'HJD[days]-2450000'),sa[2].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$')
   plt.savefig(path+'/sigma_ph.png')
   
   # [Zero_moment_WIND]
   
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
           ZM=zero_moment(x_l,y_l,H_alpha)
           zm.append(ZM)
           hjd=file.split('_')[0].replace(star,'')
           hjd=float(hjd.replace('p','.'))
           HJD_zm.append(hjd)
           if 'SONG' in file:
               ins='SONG'
           if 'MERCATOR' in file:
               ins='MERCATOR'
           if 'FIES' in file:
               ins='FIES'
           info=ZM,hjd,ins
           info=str(info)
           hdu.writelines(info)
           hdu.writelines('\n')
   X_int,Y_int=np.array(X_int),np.array(Y_int)
   X_m,Y_m=[],[]
   for j in range(len(X_int)):
       Y_m.append(Y_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
       X_m.append(X_int[j][(X_int[j]>max(M[0])) & (X_int[j]<min(M[1]))])
   X_mean,Y_mean=np.mean(X_m,0),np.mean(Y_m,0)     
   axzm.plot(X_mean,Y_mean,'k',markersize=1)
   hdu.close()
   axzm.grid()
   axzm.set_title(star+r', $H\alpha$')
   axzm.set_xlabel('Velocity [km/s]',fontsize=16),axzm.set_ylabel('Normalized flux',fontsize=16)
   axzm.axvline(x=v_alpha_1,c='g'),axzm.axvline(x=v_alpha_2,c='g')
   axzm.tick_params(axis='y', labelsize=16),axzm.tick_params(axis='x', labelsize=16)  

   plt.savefig('../stars/'+star+'/SUMMARY/'+star+'_wind.png')

   
   hdu=open(path+'/uniform_ZM.txt','r')
   ZM_old,HJD_old_zm,INS_old_zm=[],[],[]
   for j in hdu:
       dat=j.split()
       ZM_old.append(float(dat[0][1:-1]))
       HJD_old_zm.append(float(dat[1][:-1]))
       INS_old_zm.append(dat[2][:-1])
   mean_ZM=np.mean(ZM_old)
   #ZM_old=(np.array(ZM_old)-np.mean(ZM_old))/np.mean(ZM_old)
   ZM_old=(np.array(ZM_old)-np.mean(ZM_old))
   HJD_old_zm=np.array(HJD_old_zm)-2.45e6
   INS_old_zm=np.array(INS_old_zm)
   ZM_i,HJD_zm_i=[[],[],[]],[[],[],[]]
   ZM_r,HJD_zm_r=[[],[],[]],[[],[],[]]
   for zm in range(len(ZM_old)):
       if 'MERCATOR' in INS_old_zm[zm]:
           ins=0
       elif 'FIES' in INS_old_zm[zm]:
           ins=1
       elif 'SONG' in INS_old_zm[zm]:
           ins=2
       if abs(ZM_old[zm])<5*np.std(ZM_old):
           ZM_i[ins].append(ZM_old[zm]),HJD_zm_i[ins].append(HJD_old_zm[zm])
       else:
           ZM_r[ins].append(ZM_old[zm]),HJD_zm_r[ins].append(HJD_old_zm[zm])
   ZM,HJD_zm=np.concatenate(ZM_i),np.concatenate(HJD_zm_i)
   #ZM,HJD_zm,HJD_zm=ZM_old,HJD_old_zm
   np.random.seed(5)
   #pp_fm_10,sigma_fm_10=[],[]
   pp_zm_25,sigma_zm_25=[],[]
   pp_zm_50,sigma_zm_50=[],[]
   
   for i in range(25000):
       # Selection of the data to be used:
       num=[]
       rnd=np.random.randint(0,len(HJD_zm),int(0.25*len(HJD_zm)))
       for j in rnd:
           num.append(float(ZM[j]))
       #pp_fm_25.append(max(num)-min(num))
       sigma_zm_25.append(np.std(num))
       num=[]
       rnd=np.random.randint(0,len(HJD_zm),int(0.5*len(HJD_zm)))
       for j in rnd:
           num.append(float(ZM[j]))
       #pp_fm_50.append(max(num)-min(num))
       sigma_zm_50.append(np.std(num))
   sigma_zm=np.std(ZM)
   
   STD_1=[]
   step1,step10,step25,step50 = 1,2,7,50
   data_0_zm=min(HJD_zm)
   data_1_zm=max(HJD_zm)
   n_step_1,n_step_10,n_step_25=int((data_1_zm-data_0_zm)/step1)+1,int((data_1_zm-data_0_zm)/step10)+1,int((data_1_zm-data_0_zm)/step25)+1
   n_step_50=int((data_1_zm-data_0_zm)/step50)+1
   for i in range(n_step_1):
       mask1=np.logical_and(np.vstack(HJD_zm)>(i+data_0_zm),np.vstack(HJD_zm)<(i+step1+data_0_zm))
       if len(np.vstack(ZM)[mask1])>=5:
           STD_1.append(np.std(np.vstack(ZM)[mask1]))
           
   plt.close('all')
   sp,sa=plt.subplots(3)
   if len(STD_1)!=0:
       maxi,mini=max(max(STD_1),max(sigma_fm_25),max(sigma_fm_50)),min(min(STD_1),min(sigma_fm_25),min(sigma_fm_50))
   else:
       maxi,mini=max(max(sigma_fm_25),max(sigma_fm_50)),min(min(sigma_fm_25),min(sigma_fm_50))
   H_1, bins_1 = np.histogram(STD_1,bins=50,range=[mini,maxi])
   sa[0].bar((bins_1[1:]+bins_1[:-1])/2,H_1/max(H_1),alpha=0.3,width=np.mean(bins_1[1:]-bins_1[:-1]),label='1 day')
   sa[0].legend()
   H_f_25, bins_f_25 = np.histogram(sigma_zm_25,bins=50,range=[mini,maxi])
   H_m_25, bins_m_25 = np.histogram(sigma_zm_50,bins=50,range=[mini,maxi])
   sa[1].bar((bins_f_25[1:]+bins_f_25[:-1])/2,H_f_25/max(H_f_25),alpha=0.3,width=np.mean(bins_f_25[1:]-bins_f_25[:-1]),label='25%')
   sa[1].bar((bins_m_25[1:]+bins_m_25[:-1])/2,H_m_25/max(H_m_25),alpha=0.3,width=np.mean(bins_m_25[1:]-bins_m_25[:-1]),label='50%')
   sa[1].legend()
   sa[1].set_xlabel(r'$\sigma_{\langle v^0\rangle}$ [km/s]')
   sa[1].axvline(x=sigma_zm)
   sa[0].get_shared_x_axes().join(sa[0], sa[1])
   sa[2].plot(HJD_old_zm,ZM_old,'.r',markersize=1)
   sa[2].plot(HJD_zm,ZM,'.b',markersize=1)
   sa[2].grid()
   sa[2].set_xlabel(r'HJD[days]-2450000'),sa[2].set_ylabel(r'$\Delta\langle v^0\rangle$')
   plt.savefig(path+'/sigma_wnd.png')
   # [fotometria]
   if TESS==True:
       path='/home/charlie/Desktop/MSc/Grant/IAC_Science/Project/stars/'+star
       
       HJD_ph,mag=photometry(star)
       HJD_ph,mag=np.array(HJD_ph),np.array(mag)
       mag=1000*np.array(mag)
       if (max(HJD_ph)-min(HJD_ph))>30:
           #HJD_ph_1,mag_1=np.array(HJD_ph)[np.array(HJD_ph)<lim],np.array(mag)[np.array(HJD_ph)<lim]
           #HJD_ph_2,mag_2=np.array(HJD_ph)[np.array(HJD_ph)>lim],np.array(mag)[np.array(HJD_ph)>lim]
           i=1
           while min(HJD_ph[HJD_ph>(min(HJD_ph)+30*i)])-max(HJD_ph[HJD_ph<(min(HJD_ph)+30*i)])<30:
               i+=1
               if len(HJD_ph[HJD_ph>(min(HJD_ph)+30*i)])==0:
                   i-=1
                   break
           HJD_ph_1,mag_1=HJD_ph[HJD_ph<(min(HJD_ph)+30*i)],mag[HJD_ph<(min(HJD_ph)+30*i)]
           HJD_ph_2,mag_2=HJD_ph[HJD_ph>(min(HJD_ph)+30*i)],mag[HJD_ph>(min(HJD_ph)+30*i)]  
           if min(HJD_ph_2)-max(HJD_ph_1)<30:
               del(HJD_ph_2,mag_2,HJD_ph_1,mag_1)    
               HJD_ph_1,mag_1=HJD_ph,mag
       else:
           HJD_ph_1,mag_1=HJD_ph,mag
       sigma_ph_1,sigma_ph_25,sigma_ph_50=[],[],[]
       for i in range(1000):
           # Selection of the data to be used:
           num=[]
           rnd=np.random.randint(0,len(HJD_ph_1),int(0.01*len(HJD_ph_1)))
           for j in rnd:
               num.append(float(mag_1[j]))
           #pp_zm_25.append(max(num)-min(num))
           sigma_ph_1.append(np.std(num))        
           num=[]
           rnd=np.random.randint(0,len(HJD_ph_1),int(0.25*len(HJD_ph_1)))
           for j in rnd:
               num.append(float(mag_1[j]))
           #pp_zm_25.append(max(num)-min(num))
           sigma_ph_25.append(np.std(num))
           num=[]
           rnd=np.random.randint(0,len(HJD_ph_1),int(0.5*len(HJD_ph_1)))
           for j in rnd:
               num.append(float(mag_1[j]))
           sigma_ph_50.append(np.std(num))
       sigma_ph=np.std(mag)
       #plt.plot(HJD_ph_1,mag_1,'.')
       STD_1=[]
       step1= 1
       data_0_ph=min(HJD_ph_1)
       data_1_ph=max(HJD_ph_1)
       n_step_1=int((data_1_ph-data_0_ph)/step1)+1
       for i in range(n_step_1):
           mask1=np.logical_and(np.vstack(HJD_ph)>(i+HJD_ph[0]),np.vstack(HJD_ph)<(i+step1+HJD_ph[0]))
           if len(np.vstack(mag)[mask1])>=5:
               STD_1.append(np.std(np.vstack(mag)[mask1]))
       
       plt.close('all')
       sp,sa=plt.subplots(3)
       sa[2].invert_yaxis()
       maxi,mini=max(max(STD_1),max(sigma_ph_25),max(sigma_ph_50),max(sigma_ph_1)),min(min(STD_1),min(sigma_ph_25),min(sigma_ph_50),min(sigma_ph_1))
       H_1, bins_1 = np.histogram(STD_1,bins=50,range=[mini,maxi])
       sa[0].bar((bins_1[1:]+bins_1[:-1])/2,H_1/max(H_1),alpha=0.3,width=np.mean(bins_1[1:]-bins_1[:-1]),label='1 day')
       sa[0].legend()
       
       H_f_25, bins_f_25 = np.histogram(sigma_ph_25,bins=50,range=[mini,maxi])
       H_m_25, bins_m_25 = np.histogram(sigma_ph_50,bins=50,range=[mini,maxi])
       H_1, bins_1 = np.histogram(sigma_ph_1,bins=50,range=[mini,maxi])
       #sa[1].bar((bins_s_25[1:]+bins_s_25[:-1])/2,H_s_25/max(H_s_25),alpha=0.3,width=np.mean(bins_s_25[1:]-bins_s_25[:-1]),label='10%')
       sa[1].bar((bins_1[1:]+bins_1[:-1])/2,H_1/max(H_1),alpha=0.3,width=np.mean(bins_1[1:]-bins_1[:-1]),label='1% fist set')
       sa[1].bar((bins_f_25[1:]+bins_f_25[:-1])/2,H_f_25/max(H_f_25),alpha=0.3,width=np.mean(bins_f_25[1:]-bins_f_25[:-1]),label='25% fist set')
       sa[1].bar((bins_m_25[1:]+bins_m_25[:-1])/2,H_m_25/max(H_m_25),alpha=0.3,width=np.mean(bins_m_25[1:]-bins_m_25[:-1]),label='50% fist set')
       
       try:
           sigma_ph_1,sigma_ph_25,sigma_ph_50=[],[],[]
           for i in range(1000):
               # Selection of the data to be used:
               num=[]
               rnd=np.random.randint(0,len(HJD_ph_2),int(0.01*len(HJD_ph_2)))
               for j in rnd:
                   num.append(float(mag_2[j]))
               #pp_zm_25.append(max(num)-min(num))
               sigma_ph_1.append(np.std(num))        
               num=[]
               rnd=np.random.randint(0,len(HJD_ph_2),int(0.25*len(HJD_ph_2)))
               for j in rnd:
                   num.append(float(mag_2[j]))
               #pp_zm_25.append(max(num)-min(num))
               sigma_ph_25.append(np.std(num))
               num=[]
               rnd=np.random.randint(0,len(HJD_ph_2),int(0.5*len(HJD_ph_2)))
               for j in rnd:
                   num.append(float(mag_2[j]))
               sigma_ph_50.append(np.std(num))       
           H_f_25, bins_f_25 = np.histogram(sigma_ph_25,bins=50,range=[mini,maxi])
           H_m_25, bins_m_25 = np.histogram(sigma_ph_50,bins=50,range=[mini,maxi])
           H_1, bins_1 = np.histogram(sigma_ph_1,bins=50,range=[mini,maxi])
           #sa[1].bar((bins_s_25[1:]+bins_s_25[:-1])/2,H_s_25/max(H_s_25),alpha=0.3,width=np.mean(bins_s_25[1:]-bins_s_25[:-1]),label='10%')
           sa[1].bar((bins_1[1:]+bins_1[:-1])/2,H_1/max(H_1),alpha=0.3,width=np.mean(bins_1[1:]-bins_1[:-1]),label='1% second set')
           sa[1].bar((bins_f_25[1:]+bins_f_25[:-1])/2,H_f_25/max(H_f_25),alpha=0.3,width=np.mean(bins_f_25[1:]-bins_f_25[:-1]),label='25% second set')
           sa[1].bar((bins_m_25[1:]+bins_m_25[:-1])/2,H_m_25/max(H_m_25),alpha=0.3,width=np.mean(bins_m_25[1:]-bins_m_25[:-1]),label='50% second set')
       except: pass       
       sa[1].legend()
       sa[1].set_xlabel(r'$\sigma_{mag}$ [mmag]')
       sa[0].get_shared_x_axes().join(sa[0], sa[1])
       sa[1].axvline(x=sigma_ph)   
       sa[2].plot(HJD_ph,mag,'.b',markersize=1)
       sa[2].grid()
       sa[2].set_xlabel(r'HJD[days]-2450000'),sa[2].set_ylabel(r'$\Delta$Tp [mmag]')
       plt.savefig(path+'/sigma_phm.png')
       
       # Total resume 
       plt.rcParams["figure.figsize"] = (9,11)
       fig,axs=plt.subplots(3)
       
       axs[0].tick_params(axis='y', labelsize=16),axs[0].tick_params(axis='x', labelsize=16)  
       axs[1].tick_params(axis='y', labelsize=16),axs[1].tick_params(axis='x', labelsize=16)   
       axs[2].tick_params(axis='y', labelsize=16),axs[2].tick_params(axis='x', labelsize=16) 
       axs[2].invert_yaxis()
       axs[2].plot(HJD_ph,np.array(mag),'r.',markersize=1)
       axs[2].plot([], [], ' ', label=r'$\sigma_{Tp}=$'+str(round(sigma_ph,4))+' mmag')   
       # Mercator, fies, Song
       axs[0].plot(HJD_fm_r[0],FM_r[0],'rx',markersize=5)
       axs[0].plot(HJD_fm_r[1],FM_r[1],'bx',markersize=5)    
       axs[0].plot(HJD_fm_r[2],FM_r[2],'gx',markersize=5)    
       axs[0].plot(HJD_fm_i[0],FM_i[0],'r.',markersize=2)
       axs[0].plot(HJD_fm_i[1],FM_i[1],'b.',markersize=2)    
       axs[0].plot(HJD_fm_i[2],FM_i[2],'g.',markersize=2)      
       axs[0].plot([], [], ' ', label=r'$\sigma_{\langle v\rangle}=$'+str(round(sigma_fm,3))+' km/s')
       axs[1].plot(HJD_zm_r[0],ZM_r[0],'rx',markersize=5)
       axs[1].plot(HJD_zm_r[1],ZM_r[1],'bx',markersize=5)    
       axs[1].plot(HJD_zm_r[2],ZM_r[2],'gx',markersize=5)    
       axs[1].plot(HJD_zm_i[0],ZM_i[0],'r.',markersize=2)
       axs[1].plot(HJD_zm_i[1],ZM_i[1],'b.',markersize=2)    
       axs[1].plot(HJD_zm_i[2],ZM_i[2],'g.',markersize=2)   
       axs[1].plot([], [], ' ', label=r'$\sigma_{\langle v^0\rangle}=$'+str(round(sigma_zm,3))+' km/s')
       axs[0].legend(framealpha=0.1,fontsize=16),axs[1].legend(framealpha=0.1,fontsize=16),axs[2].legend(framealpha=0.1,fontsize=16)
       axs[0].axvline(x=min(HJD_ph),c='grey'),axs[0].axvline(x=max(HJD_ph),c='grey')
       axs[1].axvline(x=min(HJD_ph),c='grey'),axs[1].axvline(x=max(HJD_ph),c='grey')
       axs[0].get_shared_x_axes().join(axs[0], axs[1])
       axs[0].grid(),axs[1].grid(),axs[2].grid()
       from matplotlib.patches import Rectangle
       data_0_ph=min(HJD_ph)
       data_1_ph=max(HJD_ph)
       axs[0].get_shared_x_axes().join(axs[0], axs[1])
       axs[0].add_patch(Rectangle((data_0_fm, -sigma_fm),data_1_fm-data_0_fm, 2*sigma_fm,alpha=0.3,color='#1f77b4'))
                        #,label=r'$\sigma_{\langle v\rangle}$')
       axs[1].add_patch(Rectangle((data_0_zm, -sigma_zm),data_1_zm-data_0_zm, 2*sigma_zm,alpha=0.3,color='#1f77b4'))
                        #,label=r'$\sigma_{\langle v^0\rangle}$')
       axs[2].add_patch(Rectangle((data_0_ph, -sigma_ph),data_1_ph-data_0_ph, 2*sigma_ph,alpha=0.3,color='#1f77b4'))
                        #,label=r'$\sigma_{Tp}$')
       axs[0].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$ [km/s]',fontsize=16)    
       axs[1].set_ylabel(r'$\Delta\langle v^0\rangle$ [km/s]',fontsize=16)
       axs[2].set_ylabel(r'$\Delta$Tp [mmag]',fontsize=16)
       axs[2].set_xlabel('HJD[days]-2450000',fontsize=16)
       
       axs[0].set_title(ph_name+', '+str(line)+r'$\AA$',fontsize=16)
       axs[1].set_title(r'$H\alpha$, '+str(H_alpha)+r'$\AA$',fontsize=16)
       plt.suptitle(star,fontsize=15)
       #axs[0].legend(),axs[1].legend(),axs[2].legend()
       if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
           os.mkdir('../stars/'+star+'/SUMMARY')
       plt.savefig('../stars/'+star+'/SUMMARY/'+star+'_all.png')
   elif TESS==False:
       plt.rcParams["figure.figsize"] = (9,9)
       fig,axs=plt.subplots(2)
       axs[0].tick_params(axis='y', labelsize=16),axs[0].tick_params(axis='x', labelsize=16)  
       axs[1].tick_params(axis='y', labelsize=16),axs[1].tick_params(axis='x', labelsize=16)  
       # Mercator, fies, Song
       axs[0].plot(HJD_fm_r[0],FM_r[0],'rx',markersize=5)
       axs[0].plot(HJD_fm_r[1],FM_r[1],'bx',markersize=5)    
       axs[0].plot(HJD_fm_r[2],FM_r[2],'gx',markersize=5)    
       axs[0].plot(HJD_fm_i[0],FM_i[0],'r.',markersize=2)
       axs[0].plot(HJD_fm_i[1],FM_i[1],'b.',markersize=2)    
       axs[0].plot(HJD_fm_i[2],FM_i[2],'g.',markersize=2)      
       axs[0].plot([], [], ' ', label=r'$\sigma_{\langle v\rangle}=$'+str(round(sigma_fm,3))+' km/s')
       axs[1].plot(HJD_zm_r[0],ZM_r[0],'rx',markersize=5)
       axs[1].plot(HJD_zm_r[1],ZM_r[1],'bx',markersize=5)    
       axs[1].plot(HJD_zm_r[2],ZM_r[2],'gx',markersize=5)    
       axs[1].plot(HJD_zm_i[0],ZM_i[0],'r.',markersize=2)
       axs[1].plot(HJD_zm_i[1],ZM_i[1],'b.',markersize=2)    
       axs[1].plot(HJD_zm_i[2],ZM_i[2],'g.',markersize=2)   
       axs[1].plot([], [], ' ', label=r'$\sigma_{\langle v^0\rangle}=$'+str(round(sigma_zm,3)))
       axs[0].legend(framealpha=0.1),axs[1].legend(framealpha=0.1)
       axs[0].get_shared_x_axes().join(axs[0], axs[1])
       axs[0].grid(),axs[1].grid()
       from matplotlib.patches import Rectangle
       axs[0].get_shared_x_axes().join(axs[0], axs[1])
       axs[0].add_patch(Rectangle((data_0_fm, -sigma_fm),data_1_fm-data_0_fm, 2*sigma_fm,alpha=0.3,color='#1f77b4'))
                        #,label=r'$\sigma_{\langle v\rangle}$')
       axs[1].add_patch(Rectangle((data_0_zm, -sigma_zm),data_1_zm-data_0_zm, 2*sigma_zm,alpha=0.3,color='#1f77b4'))
                        #,label=r'$\sigma_{\langle v^0\rangle}$')
                        #,label=r'$\sigma_{Tp}$')
       axs[0].set_ylabel(r'$\langle v\rangle-\langle v\rangle_0$ [km/s]')    
       axs[1].set_ylabel(r'$\Delta\langle v^0\rangle$ [km/s]')
       axs[0].set_title(ph_name+', '+str(line)+r'$\AA$')
       axs[1].set_title(r'$H\alpha$, '+str(H_alpha)+r'$\AA$')
       axs[1].set_xlabel('HJD[days]-2450000')
       plt.suptitle(star,fontsize=15)
       #axs[0].legend(),axs[1].legend(),axs[2].legend()
       if os.path.exists('../stars/'+star+'/SUMMARY')!=True:
           os.mkdir('../stars/'+star+'/SUMMARY')
       plt.savefig('../stars/'+star+'/SUMMARY/'+star+'_all.png')
       sigma_ph=np.nan
   return star,sigma_fm,sigma_zm,mean_ZM,sigma_ph,len(FM),len(ZM)
 
   