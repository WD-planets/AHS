'''
This code has different functions in terms of how many dips you see per transit. 
It's rudimentary. I tried to generalise it with lines 38 to 54 but failed. If you want
to explore that option what fails is the tuple to feed to curve fit

Procedure: I let the out of transit level, the depth factor scale, the
time at minimum flux and the ingress and egress time-scales as free parameters for the first time we fit 
a transit. If I want to track down a transit in different orbits: I fix the depth factor, ingress and 
egress and just let the out-of transit and time at minimum flux as free parameters 
(in order to "preserve") the shape of the transit.

The numbers in the first brackets are the values, and the second bracket correspond 
to uncertainties of the fit.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
from PyAstronomy.pyasl import binningx0dt

#This is the single and simplest case AHS with one dip
def ahs_s(t,oof,a,t0,ti,te):
    # oof =  out of transit flux
    # a = depth factor scale
    # t0 = time at minimum flux
    # ti, te = ingress and egress timescales 
    return(oof/((a/(np.exp((t-t0)/te)+np.exp((t0-t)/ti))+1.0)))

#two,three and four dips inside same transit
def ahs2(t,oof,a1,t01,ti1,te1,a2,t02,ti2,te2):
    return( oof/((a1/(np.exp((t-t01)/te1)+np.exp((t01-t)/ti1))+1.0)*(a2/(np.exp((t-t02)/te2)+np.exp((t02-t)/ti2))+1.0)) )    

def ahs3(t,oof,a1,t01,ti1,te1,a2,t02,ti2,te2,a3,t03,ti3,te3):
    return( oof/((a1/(np.exp((t-t01)/te1)+np.exp((t01-t)/ti1))+1.0)*(a2/(np.exp((t-t02)/te2)+np.exp((t02-t)/ti2))+1.0)*(a3/(np.exp((t-t03)/te3)+np.exp((t03-t)/ti3))+1.0)) )    

def ahs4(t,oof,a1,t01,ti1,te1,a2,t02,ti2,te2,a3,t03,ti3,te3,a4,t04,ti4,te4):
    return( oof/((a1/(np.exp((t-t01)/te1)+np.exp((t01-t)/ti1))+1.0)*(a2/(np.exp((t-t02)/te2)+np.exp((t02-t)/ti2))+1.0)*(a3/(np.exp((t-t03)/te3)+np.exp((t03-t)/ti3))+1.0)*(a4/(np.exp((t-t04)/te4)+np.exp((t04-t)/ti4))+1.0)) )    

#The following two functions were to try to generalise things, but I can't make them work
#This is just the "productory"
def product_recursive(numbers):
    if not numbers:
        return 1 # base case: list is empty
    return numbers[0] * product_recursive(numbers[1:]) # recursive case: multiply first element by product of the rest of the list

#This is the function I failed to make work because tuples are not allowed to be fed to curve_fit
def ahs_gral(t,oof,a,t0,ti,te):
    n = len(a)
    each = []
    for i in range(n):
        print('again',a[i],t0[i],ti[i],te[i])
        each.append(a[i]/(np.exp((t-t0[i])/te[i])+np.exp((t0[i]-t)/ti[i]))+1.0 )
        print(each)
    print(product_recursive(each))
    return(oof/product_recursive(each))


#I don't remember what these things do (they are related to the specific problem I was solving at that time)
def conver(zu):
    return((zu-0.39)*24*60-12.8664)
    
def tomin(q):
    return(q*24*60)

day1=np.loadtxt('/Users/paulaizquierdo/Documents/Astrofisica/DOC/WD1054/ucam_mar28.dat')

bjd    = day1[:,0] - day1[0,0]
star   = day1[:,2]
comp   = day1[:,6]
erstar = day1[:,4]
ercomp = day1[:,8]

day00 = 2458567.49298

flux   = star/comp
erflux = np.sqrt( ((ercomp*star/(comp**2))**2) + ((erstar/comp)**2) ) 

filt_med = (bjd>0.0246) & (bjd<0.029)
data_bin, peo = binningx0dt(bjd,flux/np.mean(flux[filt_med]),yerr=erflux,reduceBy=2, useBinCenter=True)

#filtering the data belonging to the transit
transit = (data_bin[:,0] > 0.0962) & (data_bin[:,0] < 0.1043)

#plotting the whole LC
plt.figure()
plt.plot(data_bin[:,0],data_bin[:,1],'-')
plt.xlim(0,bjd[-1])
plt.xlabel('BJD')
plt.ylabel('Flux')
plt.show()

# curve_fit(function, x_array, y_array, initial guesses, sigma)
aj,cov = sc.curve_fit(ahs2,data_bin[transit,0],data_bin[transit,1],(0.97,0.05,0.097,1e-4,1e-3,0.04,0.102,1e-4,1e-4),sigma=data_bin[transit,2])
sd = np.sqrt(np.diag(cov))

f=data_bin[transit,0]*0.+1.0 #Array of ones with the dimension of the transit
fit = f*ahs2(data_bin[transit,0],*aj) #best fit to plot

print('parameters and uncertaintues i=2',aj,sd)
print(aj[2]+(day1[0,0]-day00),aj[6]+(day1[0,0]-day00))


plt.figure(2)
plt.plot(data_bin[transit,0],data_bin[transit,1],'k-') 
plt.plot(data_bin[transit,0],data_bin[transit,1],'ko') 
plt.errorbar(data_bin[transit,0],data_bin[transit,1],yerr=data_bin[transit,2],fmt='', ecolor=None)
plt.plot(data_bin[transit,0],fit,'r-')
plt.xlabel('BJD')
plt.ylabel('Flux')
plt.show()

