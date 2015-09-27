#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import scipy as sci
import math as m
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib as mpl

#Function to calculate refractive indices of Gold
def nDRUDE(freq, wp, gamma):
    damping=np.vectorize(complex)(0,freq*gamma)
    return np.sqrt(1-((wp**2)/(freq**2+damping)))

#Function to calculate Fresnel reflection coefficients for both interfaces, in the s and p planes
def fresnelr(q,r1,r2,r3,wav,l,pl,interface):
    znondim1=np.sqrt(r1**2-q**2)
    znondim2=np.sqrt(r2**2-q**2)
    znondim3=np.sqrt(r3**2-q**2)
    
    if interface=='onetwo':
        if pl=='s':
            return (znondim1-znondim2)/(znondim1+znondim2)
        elif pl=='p':
            return -((znondim2/(r2**2))-(znondim1/(r1**2)))/((znondim2/(r2**2))+(znondim1/(r1**2)))
    elif interface=='twothree':
        if pl=='s':
            return (znondim2-znondim3)/(znondim2+znondim3)
        elif pl=='p':
            return -((znondim3/(r3**2))-(znondim2/(r2**2)))/((znondim3/(r3**2))+(znondim2/(r2**2)))
        
#Function to calculate composite reflection coefficients in the s and p planes from fresnel reflection coefficients
def compr(q,r1,r2,r3,wav,l,pla):
    pow=np.vectorize(complex)(0, 2*2*m.pi*(1/wav)*np.sqrt(r2**2-q**2)*l)
    refls12=fresnelr(q,r1,r2,r3,wav,l,'s','onetwo')
    reflp12=fresnelr(q,r1,r2,r3,wav,l,'p','onetwo')
    refls23=fresnelr(q,r1,r2,r3,wav,l,'s','twothree')
    reflp23=fresnelr(q,r1,r2,r3,wav,l,'p','twothree')
    
    if (pla=='s'):
        return (refls12+refls23*np.exp(pow))/(1+refls12*refls23*np.exp(pow))
    elif (pla=='p'):
        return (reflp12+reflp23*np.exp(pow))/(1+reflp12*reflp23*np.exp(pow))

#Function to calculate angular emissivity, to be included in the next function for integration
def angularem(q,r1,r2,r3,wav,l):
    reflscomp=compr(q,r1,r2,r3,wav, l,'s')
    reflpcomp=compr(q,r1,r2,r3,wav, l,'p')
    return (1-(np.abs(reflscomp))**2)+(1-(np.abs(reflpcomp))**2)
    
#Function to return the expression to integrate across all angles 
def integrand(q,r1,r2,r3,wav,l):
    return q*angularem(q,r1,r2,r3,wav,l)

#Integrate angular emissivity across all anglues to get hemispherical emissivity 
def spechemem(r1,r2,r3,wav,l):
    she=quad(integrand,0,1,args=(r1,r2,r3,wav,l))
    return she[0]

#Function to be called from MAIN file
def generate_sel_em_emissivity(filename):
    
    if(filename=="SiO2data.txt"):
        material="SiO2-Au"
    elif(filename=="BNdata.txt"):
        material="BN-Au"    
    
    #Refractive Index of Vacuum
    n1=1.0
    
    #Read refractive indices for SiO2
    n2values=np.genfromtxt(filename, dtype=float)
    
    #Read n and k values from array "n2values" and store into an array in n+ik form
    n2old=np.vectorize(complex)(n2values[:,1],1*n2values[:,2])
    
    #Create function f to interpolate for n2 values
    f=interpolate.interp1d(n2values[:,0],n2old)
    
    #Define wavelength intervals
    wavelengths=np.arange(2.6,32.1,0.1,dtype=float)
    
    #Get n2 values for these intervals by interpolating
    n2=f(wavelengths)
    
    #Define constants to calculate frequency in eV
    echarge = 1.6e-19
    planck = 6.626e-34
    cvac = 3e8
    
    #Calculate frequencies in eV
    wavelengthinmicrons = wavelengths
    freqinev = planck*cvac/(echarge*wavelengthinmicrons*1e-6)
        
    #Calculate refractive indices for Gold
    n3=nDRUDE(freqinev,9.03,0.0267)
    
    #Create array to store hemispherical emissivity values    
    hem=np.zeros(shape=(len(n2),4),dtype=float)
    
    #Store emissivity values into array
    for x in range(0,len(n2)):
        l1=0.2
        l2=0.5
        l3=1.0
        r1=n1
        r2=n2[x]
        r3=n3[x]
        wav=wavelengths[x]    
        hem[x,0]=wavelengths[x]
        hem[x,1]=spechemem(r1,r2,r3,wav,l1)
        hem[x,2]=spechemem(r1,r2,r3,wav,l2)
        hem[x,3]=spechemem(r1,r2,r3,wav,l3)
        print hem[x,:]
    
    #Save hemispherical emissivity data to text file for use in Problem 3
    np.savetxt("Selem_Hem_Emissivityfor%s.txt" % material,hem,delimiter=' ')
    
    
    #Plot
    fig=plt.figure()
    fig.suptitle("Hem. Emissivities of %s Selective Emitter" % material)
    ax1=fig.add_subplot(221)
    ax1.set_ylabel("Emissivity (R-200nm|B-500nm|G-1000nm)")
    ax1.set_xlabel("Wavelength (microns)")
    ax1.set_xlim([2,33]) 
    ax1.set_ylim([0,1])
    ax1.set_title("%s (all thicknesses)" % material)
    ax1.plot(hem[:,0],hem[:,1],color='r',label="Emissivity vs. Wavelength",linewidth=2)
    ax5=plt.twinx(ax1)
    ax5.set_ylim([0,1])
    ax5.set_xlim([2,33]) 
    ax5.plot(hem[:,0],hem[:,2],color='b',linewidth=2)
    ax6=plt.twinx(ax1)
    ax6.set_ylim([0,1])
    ax6.set_xlim([2,33]) 
    ax6.plot(hem[:,0],hem[:,3],color='g',linewidth=2)
    
    ax2=fig.add_subplot(222)
    ax2.set_ylabel("Emissivity (200nm)")
    ax2.set_xlabel("Wavelength (microns)")
    ax2.set_xlim([2,33]) 
    ax2.set_ylim([0,1])
    ax2.set_title("%s (200nm)" % material)
    ax2.plot(hem[:,0],hem[:,1],color='r',label="Emissivity vs. Wavelength",linewidth=2)
 
    ax3=fig.add_subplot(223)
    ax3.set_ylabel("Emissivity (500nm)")
    ax3.set_xlabel("Wavelength (microns)")
    ax3.set_xlim([2,33]) 
    ax3.set_ylim([0,1])
    ax3.set_title("%s (500nm)" % material)
    ax3.plot(hem[:,0],hem[:,2],color='b',label="Emissivity vs. Wavelength",linewidth=2)
    
    ax4=fig.add_subplot(224)
    ax4.set_ylabel("Emissivity (1000nm)")
    ax4.set_xlabel("Wavelength (microns)")
    ax4.set_xlim([2,33]) 
    ax4.set_ylim([0,1])
    ax4.set_title("%s (1000nm)" % material)
    ax4.plot(hem[:,0],hem[:,3],color='g',label="Emissivity vs. Wavelength",linewidth=2)
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig("Emissivityfor%s.png" % material)
