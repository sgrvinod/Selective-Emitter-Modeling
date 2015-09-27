#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import scipy as sci
import math as m
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate
import os

#Function to calculate Fresnel reflection coefficients in the s and p planes
def fresnelr(q,r1,r2,pl):
    znondim1=np.sqrt(r1**2-q**2)
    znondim2=np.sqrt(r2**2-q**2)
            
    if pl=='s':
        return (znondim1-znondim2)/(znondim1+znondim2)
    elif pl=='p':
        return -((znondim2/(r2**2))-(znondim1/(r1**2)))/((znondim2/(r2**2))+(znondim1/(r1**2)))
   
#Function to calculate angular emissivity, to be included in the next function for integration
def angularem(q,r1,r2):
    refls=fresnelr(q,r1,r2,'s')
    reflp=fresnelr(q,r1,r2,'p')
    return (1-(np.abs(refls))**2)+(1-(np.abs(reflp))**2)
    
#Function to return the expression to be integrated across all angles 
def integrand(q,r1,r2):
    return q*angularem(q,r1,r2)

#Integrate angular emissivities across all anglues to get hemispherical emissivity 
def spechemem(r1,r2):
    she=quad(integrand,0,1,args=(r1,r2))
    return she[0]

#Fuction to be called in MAIN file
def generate_singlelayerem(filename):
    
    if(filename=="SiO2data.txt"):
        material="SiO2"
    elif(filename=="BNdata.txt"):
        material="BN"    
    
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
    
    #Create array to store hemispherical emissivity values    
    hem=np.zeros(shape=(len(n2),2),dtype=float)
    
    #Store hemispherical emissivity values into array
    for x in range(0,len(n2)):
        r1=n1
        r2=n2[x]
        hem[x,0]=wavelengths[x]
        hem[x,1]=spechemem(r1,r2)
        print hem[x,:]
    
    #Save hemispherical emissivity data to text file for use in Problem 3
    np.savetxt("Singlelayer_Hem_Emissivityfor%s.txt" % material,hem,delimiter=' ')
    
    #Plot
    fig=plt.figure()
    fig.suptitle("Hem. Emissivities of %s" % material)
    ax1=fig.add_subplot(111)
    ax1.set_ylabel("Emissivity")
    ax1.set_xlabel("Wavelength (microns)")
    ax1.set_xlim([2,33])    
    ax1.set_ylim([0,1])
    ax1.plot(hem[:,0],hem[:,1],color='k',label="Emissivity vs. Wavelength",linewidth=2)
    fig.savefig("Hem_Emissivityfor%s.png" % material)


