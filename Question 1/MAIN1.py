#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import scipy as sci
import math as mat
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib
import matplotlib.pyplot as plt

#Linear equation definition
def linearfit(x, d, e):
    return d*x+e

#Quadratic equation definition
def quadraticfit(x, d, e, n):
    return ((d*(x**2)+e*x+n))

#Integrand for calculating hemispherical values of transmissivity
def integrandforhem(x):
    return x

#Function to model atmospheric transmissivities at specified humidity
def generatehemtransdata(filename):
    
    if(filename=="2.3mm.txt"):
        humidity="2,3mm"
    elif(filename=="4.3mm.txt"):
        humidity="4,3mm"
    elif(filename=="7.6mm.txt"):
        humidity="7,6mm"
    elif(filename=="10mm.txt"):
        humidity="10,00mm"
    
    #Read transmissivity values for 2.3mm from file into array "trans"
    trans=np.genfromtxt(filename, dtype=float)
    
    #Calculate length of required simplified array
    lensimpl=(320-26)+1
    
    #Create this simplified array "transsimpl"
    transsimpl=np.zeros(shape=(lensimpl,2),dtype=float)
    
    #Create a count for keeping track of indices of "transsimpl"
    a=0
    
    #For convenience, express wavelength in 10^-1 microns instead of microns
    #note:We are doing this because we want to extract wavelength multiples of 0.1 using the % operator
    for i in range(len(trans)):
        trans[i,0]=trans[i,0]*10
    
    #Store transmissivity values between 2.6 and 5.5 microns into "transsimpl" from "trans"
    for i in range(len(trans)):
        if (trans[i,0]%1==0):
            if(trans[i,0]<26):
                continue
            if(trans[i,0]>25 and trans[i,0]<56):
                transsimpl[a,0]=trans[i,0]
                transsimpl[a,1]=trans[i,1]
                a=a+1
                
    #Store transmissivity=0 between 5.6 and 6.9 microns into "transsimpl"
    w=56
    for i in range(a,a+14):
        transsimpl[i,0]=w
        transsimpl[i,1]=0
        w=w+1
    
    #Store transmissivity values between 7.0 and 26.0 microns into "transsimpl" from "trans"
    a=a+14
    for i in range(len(trans)):
        if (trans[i,0]%1==0):
            if(trans[i,0]>260):
                continue
            if(trans[i,0]>69 and trans[i,0]<261):
                transsimpl[a,0]=trans[i,0]
                transsimpl[a,1]=trans[i,1]
                a=a+1
    
    #Store transmissivity=0.1 between 26.1 and 32.0 microns into "transsimpl"
    w=261
    for i in range(a,a+60):
        transsimpl[i,0]=w
        transsimpl[i,1]=0.1
        w=w+1
    
    #Convert wavelengths from 10^-1 microns back to microns
    for i in range(len(transsimpl)):
        transsimpl[i,0]=transsimpl[i,0]/10
        
    #Create a new array "Transapprox" to store approximated transmissivity values
    transapprox=np.zeros(shape=(lensimpl,2),dtype=float)
    
    #Calculate indices for intervals that must be curve fit
    for i in range(len(transsimpl)):
        if transsimpl[i,0]==3.7:
            int1st=i
        elif transsimpl[i,0]==4.0:
            int1stop=i
            int2st=i
        elif transsimpl[i,0]==4.2:
            int2stop=i
        elif transsimpl[i,0]==4.4:
            int3st=i
        elif transsimpl[i,0]==4.6:
            int3stop=i
        elif transsimpl[i,0]==4.9:
            int4st=i
        elif transsimpl[i,0]==5.6:
            int4stop=i
        elif transsimpl[i,0]==7.0:
            int5st=i
        elif transsimpl[i,0]==8.0:
            int5stop=i
        elif transsimpl[i,0]==13.0:
            int6st=i
        elif transsimpl[i,0]==14.2:
            int6stop=i
        elif transsimpl[i,0]==15.8:
            int7st=i
        elif transsimpl[i,0]==18.0:
            int7stop=i
            int8st=i
        elif transsimpl[i,0]==26.0:
            int8stop=i
            
    #Get parameters of linear equation defined above - d, e - for required interval
    def calculatelinearfit(f,g):
        x=transsimpl[f:g,0]
        y=transsimpl[f:g,1]
        popt,pcov=curve_fit(linearfit,x,y)
        return popt
    
    #Get parameters of quadratic equation defined above - d, e, n - for required interval
    def calculatequadraticfit(f,g):
        x=transsimpl[f:g,0]
        y=transsimpl[f:g,1]
        popt,pcov=curve_fit(quadraticfit,x,y)
        return popt
        
    #Get box value for required interval
    def boxvalue(b,c):
        for i in range(len(transsimpl)):
            if transsimpl[i,0]==b:
                st=i
            if transsimpl[i,0]==c:
                stop=i
        return (np.trapz(transsimpl[st:stop,1],dx=0.1))/(c-b)
    
    #Fit to curves using curve parameters or box values
    for i in range(len(transsimpl)):
        #For 2.6µm-2.9µm, T=0
        if (transsimpl[i,0]>=2.6 and transsimpl[i,0]<=2.9):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=0
        #For 2.9µm-3.7µm, T=boxvalue
        elif (transsimpl[i,0]>=2.9 and transsimpl[i,0]<=3.7):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=boxvalue(2.9,3.7)
        #For 3.7µm-4.0µm, T is fitted linearly
        elif (transsimpl[i,0]>=3.7 and transsimpl[i,0]<=4.0):
            transapprox[i,0]=transsimpl[i,0]
            h,j=calculatelinearfit(int1st+1,int1stop+1)
            transapprox[i,1]=h*transapprox[i,0]+j
        #For 4.0µm-4.2µm, T is fitted linearly
        elif (transsimpl[i,0]>=4.0 and transsimpl[i,0]<=4.2):
            transapprox[i,0]=transsimpl[i,0]
            k,l=calculatelinearfit(int2st+1,int2stop+1)
            transapprox[i,1]=k*(transapprox[i,0])+l
        #For 4.2µm-4.4µm, T=0
        elif (transsimpl[i,0]>=4.2 and transsimpl[i,0]<=4.4):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=0
        #For 4.4µm-4.6µm, T is fitted linearly
        elif (transsimpl[i,0]>=4.4 and transsimpl[i,0]<=4.6):
            transapprox[i,0]=transsimpl[i,0]
            m,p=calculatelinearfit(int3st+1,int3stop+1)
            transapprox[i,1]=m*(transapprox[i,0])+p
        #For 4.6µm-4.9µm, T=boxvalue
        elif (transsimpl[i,0]>=4.6 and transsimpl[i,0]<=4.9):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=boxvalue(4.6,4.9)
        #For 4.9µm-5.6µm, T is fitted linearly
        elif (transsimpl[i,0]>=4.9 and transsimpl[i,0]<=5.6):
            transapprox[i,0]=transsimpl[i,0]
            q,r=calculatelinearfit(int4st+1,int4stop+1)
            transapprox[i,1]=q*(transapprox[i,0])+r
        #For 5.6µm-7.0µm, T=0
        elif (transsimpl[i,0]>=5.6 and transsimpl[i,0]<=7.0):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=0
        #For 7.0µm-8.0µm, T is fitted linearly
        elif (transsimpl[i,0]>=7.0 and transsimpl[i,0]<=8.0):
            transapprox[i,0]=transsimpl[i,0]
            s,t=calculatelinearfit(int5st+1,int5stop+1)
            transapprox[i,1]=s*(transapprox[i,0])+t
        #For 8.0µm-13.0µm, T=boxvalue
        elif (transsimpl[i,0]>=8.0 and transsimpl[i,0]<=13.0):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=boxvalue(8.0,13.0)
        #For 13.0µm-14.2µm, T is fitted linearly
        elif (transsimpl[i,0]>=13.0 and transsimpl[i,0]<=14.2):
            transapprox[i,0]=transsimpl[i,0]
            u,v=calculatelinearfit(int6st+1,int6stop+1)
            transapprox[i,1]=u*(transapprox[i,0])+v
        #For 14.2µm-15.8µm, T=0
        elif (transsimpl[i,0]>=14.2 and transsimpl[i,0]<=15.8):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=0
        #For 15.8µm-18.0µm, T is fitted to a quadratic
        elif (transsimpl[i,0]>=15.8 and transsimpl[i,0]<=18.0):
            transapprox[i,0]=transsimpl[i,0]
            q1,q2,q3=calculatequadraticfit(int7st+1,int7stop+1)
            transapprox[i,1]=q1*(transapprox[i,0]**2)+q2*transapprox[i,0]+q3
        #For 18.0µm-26.0.0µm, T is fitted linearly
        elif (transsimpl[i,0]>=18.0 and transsimpl[i,0]<=26.0):
            transapprox[i,0]=transsimpl[i,0]
            l1,l2=calculatelinearfit(int8st+1,int8stop+1)
            transapprox[i,1]=l1*(transapprox[i,0])+l2
        #For 26.0µm-32.0µm, T=0.1
        elif (transsimpl[i,0]>=26.0 and transsimpl[i,0]<=32.0):
            transapprox[i,0]=transsimpl[i,0]
            transapprox[i,1]=0.1
    
    #Create new arrays to hold spectral hemispherical transmissivity values
    transspecsimpl=np.zeros(shape=(lensimpl,2),dtype=float)
    transspecapprox=np.zeros(shape=(lensimpl,2),dtype=float)
    
    #Integrate integrand above
    hemfactor=quad(integrandforhem,0,1)
    
    #Compute hemispherical values of transmissivity
    for i in range(len(transapprox)):
        transspecsimpl[i,0]=transsimpl[i,0]
        transspecapprox[i,0]=transsimpl[i,0]
        
        transspecsimpl[i,1]=2*transsimpl[i,1]*hemfactor[0]
        transspecapprox[i,1]=2*transapprox[i,1]*hemfactor[0]
        
    #Save Exact Hemispherical Transmissivity Data to text file for use in Question 3
    np.savetxt("Exact_Hem_Transfor%s.txt" % humidity,transspecsimpl,delimiter=' ')
    
    #Save Approximated HemisphericalT ransmissivity Data to text file for use Question 3
    np.savetxt("Approx_Hem_Transfor%s.txt" % humidity,transspecapprox,delimiter=' ')
    
    #Plot approximated hemispherical transmissivity values against wavelength
    fig=plt.figure()
    ax1=fig.add_subplot(111)    
    ax1.set_label(humidity)
    ax1.set_ylabel("Transmissivity - Approx (Red)")
    ax1.set_xlabel("Wavelength (microns)")
    ax1.set_ylim([0,1])
    ax1.set_xlim([2,32])
    ax1.set_title("Hemispherical Transmissivity - %s" % humidity)
    ax1.plot(transspecapprox[:,0],transspecapprox[:,1],color='r',label="Approx. Transmissivity vs. Wavelength",linewidth=2)
    
    #Overlay exact transmissivity values
    ax2=plt.twinx(ax1)
    ax2.set_ylabel('Transmissivity - Exact (Blue)')
    ax2.set_ylim([0,1])
    ax2.plot(transspecsimpl[:,0],transspecsimpl[:,1],color='b',linewidth=2)
    fig.savefig("Hem_Emissivityfor%s.png" % humidity)
    
#generatehemtransdata has been defined
    
generatehemtransdata("2.3mm.txt")

generatehemtransdata("4.3mm.txt")

generatehemtransdata("7.6mm.txt")

generatehemtransdata("10mm.txt")

    
    
    
         
        
        
        
        
        
            
        
    
    
            
                
       
    
    
    
            
    
        
            
            
        
       
