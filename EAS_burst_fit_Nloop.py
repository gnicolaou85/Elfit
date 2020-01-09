import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from pylab import *
from optparse import OptionParser
from matplotlib import ticker
import scipy
from scipy import optimize
import time

### Constants #####
me=9.11*10.**(-31.)
joule=1.602*10.**(-19.)
rad=(np.pi/180.)
cmap='Greens'
###################

####   INTRUMENT's PARAMETERS #########################################
elevation_table_min=np.arange(16)*6.-45.0
elevation_table_max=elevation_table_min+6.
elevation_table_middle=0.5*(elevation_table_min+elevation_table_max)

elevation_plot=np.zeros(17)
elevation_plot[0:16]=elevation_table_min[:]
elevation_plot[-1]=elevation_table_max[-1]

energy_table_middle=10.**(np.linspace(np.log10(1),np.log10(5000),64))
energy_table_min=energy_table_middle-0.0625*energy_table_middle
energy_table_max=energy_table_middle+0.0625*energy_table_middle                                         
                                                      
energy_plot=np.zeros(65)
energy_plot[0:64]=energy_table_min[:]
energy_plot[-1]=energy_table_max[-1]

azimuth_table_min=np.arange(32)*11.25
azimuth_table_max=np.arange(32)*11.25+11.25
azimuth_table_middle=0.5*(azimuth_table_min+azimuth_table_max)

azimuth_plot=np.zeros(33)
azimuth_plot[0:32]=azimuth_table_min[:]
azimuth_plot[-1]=azimuth_table_max[-1]

Geometric_factor=( 1.6*10.**(-8.) ) * (0.5) * (11.25/15.)  ## we scale LEEA G factor by 2 and then adjust to azimuth resolution
DT=0.00087   ## acquisition time


                             ####   grid ####
elevation,energy,azimuth=np.meshgrid(elevation_table_middle,energy_table_middle,azimuth_table_middle)

U=np.sqrt(2.*energy*joule/me)
Ux=U*np.cos(elevation*rad)*np.cos(azimuth*rad)
Uy=U*np.cos(elevation*rad)*np.sin(azimuth*rad)
Uz=U*np.sin(elevation*rad)
#######################################################################################################

######################
## Input Parameters ##
######################                                                                           ##########################################################################################
Density_array=(10.**linspace(np.log10(5),np.log10(500),20)) * 1000000. # density in cm^-3        ## modify here the range of density to examine                                          ##
                                                                                                 ## e.g., in Nicolaou et al. we investiogated n=7cm-3 and a range from 5cm-3 to 500cm-3  ##
#print (Density_array)                                                                           ##########################################################################################
#quit()

E0=0.7115

anisotropy=2.
Eth_par=10.
Eth_perp=anisotropy*Eth_par

uth_par=np.sqrt(2.*Eth_par*joule/me)
uth_perp=np.sqrt(2.*Eth_perp*joule/me)

kappa=3.
alpha=45.     # angle between velocity (x-axis) and B_field (on top hat plane)


###############################
###############################
U0=np.sqrt(2.*E0*joule/me)
Ux0=U0                       ## Bulk velocity fixed along x-axis same for B-field  
Uy0=0.
Uz0=0.

U0_par=Ux0*np.cos(alpha*rad)+Uy0*np.sin(alpha*rad)
U0_perp=Ux0*np.sin(alpha*rad)-Uy0*np.cos(alpha*rad)         ### U perp on top hat plane. (Uz is also Uperp)
####################################
####################################
kappa_results_array=np.zeros([len(Density_array),200])
Density_results_array=np.zeros([len(Density_array),200])
Tpar_results_array=np.zeros([len(Density_array),200])
Tperp_results_array=np.zeros([len(Density_array),200])
exp_counts_array=np.zeros(len(Density_array))
for Density_deik in range(len(Density_array)):
  Density=Density_array[Density_deik]
  print (Density)
  for sample_deik in range(200):

    ## building the INPUT distribution ##
    Upar=Ux*np.cos(alpha*rad)+Uy*np.sin(alpha*rad)
    Uperp= Ux*np.sin(alpha*rad)-Uy*np.cos(alpha*rad)  ### U perp on top hat plane. (Uz is also Uperp)

    parenthesis=(Upar-U0_par)**2./(uth_par**2.) + (Uperp-U0_perp)**2./(uth_perp**2.)+ (Uz**2.)/(uth_perp**2.)
          
    norm=np.pi**(-1.5)*(uth_par**(-1.))*(uth_perp**(-2.))*((kappa-1.5)**(-1.5))*(math.gamma(kappa+1.)/math.gamma(kappa-0.5))
    f=Density*norm*((1.+parenthesis/(kappa-1.5))**(-kappa-1.))

    C=f*(U**4.)*Geometric_factor*DT   ## counts
    exp_counts_array[Density_deik]=np.max(C)  
    C_out=np.random.poisson(C)                    #actual_measurement in 3D
    C_out_burst=C_out[:,7,:]
    #print (elevation_table_middle[7])
    #quit()

##    ###############################################
##    ## plot the counts in the instrument's frame ##
##    ###############################################
##    fig=plt.figure(figsize=(6,6))
##    ax1 = fig.add_axes([0.15,0.25,0.6,0.6], projection='polar')
##    ax1.set_rticks([1, 2, 3, 4])
##    ax1.tick_params(labelsize=12)
##    spec=ax1.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(C_out_burst),cmap='Greens')
##    cax1=fig.add_axes([0.15,0.13,0.6,0.04])
##    cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(Counts)',cax=cax1)#,pad=0.08,fraction=0.054)
##    ax1.set_xlim([0.*rad,360.*rad])
##    ax1.grid(linestyle='--',alpha=0.5)
##    plt.show()
##    ###############################################
##    ###############################################

    ###########################
    ## getting f from counts ##
    ###########################
    f_out=C_out/(Geometric_factor*(U**(4.))*DT)          ###########
    f_out_burst=f_out[:,7,:]

    Ux_burst=Ux[:,7,:]
    Uy_burst=Uy[:,7,:]
    Uz_burst=Uz[:,7,:]
    ###########################
    ###########################

##    #####################################################
##    ## plot the distribution in the instrument's frame ##
##    #####################################################
##    fig=plt.figure(figsize=(6,6))
##    ax1 = fig.add_axes([0.15,0.25,0.6,0.6], projection='polar')
##    ax1.set_rticks([1, 2, 3, 4])
##    ax1.tick_params(labelsize=12)
##    spec=ax1.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(f_out_burst),cmap='Greens')
##    cax1=fig.add_axes([0.15,0.13,0.6,0.04])
##    cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(f_out)',cax=cax1)#,pad=0.08,fraction=0.054)
##    ax1.set_xlim([0.*rad,360.*rad])
##    ax1.grid(linestyle='--',alpha=0.5)
##    plt.show()
##    #####################################################
##    #####################################################

    #############
    ## Fitting ##
    #############
    def fopt(p):
          
      Nf=p[0]
      uth_parf=p[1]
      uth_perpf=p[2]
      kappaf=p[3]
      
      parf=(Upar-U0_par)**2./(uth_parf**2.) + (Uperp-U0_perp)**2./(uth_perpf**2.)+ (Uz**2.)/(uth_perpf**2.)
              
      normf=np.pi**(-1.5)*(uth_parf**(-1.))*(uth_perpf**(-2.))*((kappaf-1.5)**(-1.5))*(math.gamma(kappaf+1.)/math.gamma(kappaf-0.5))
      f_fit=Nf*normf*((1.+parf/(kappaf-1.5))**(-kappaf-1.))

      C_fit=f_fit*(U**4.)*Geometric_factor*DT

      C_fit_burst=C_fit[:,7,:]

      #plt.subplot(1,2,1)
      #sp1=plt.pcolor(C_out_burst)
      #cb1=plt.colorbar(sp1)
      #plt.subplot(1,2,2)
      #sp2=plt.pcolor(C_fit_burst)
      #cb2=plt.colorbar(sp2)
      #plt.show()                                                     #############################################################
      variance=1.*C_out_burst                                         ## comment and modify accordingly for fitting C<1 and C>=1 ##
      deik=np.where(C_out_burst>=1.)                                   #############################################################
      #number_of_points=len(C_out_burst)*len(C_out_burst[0])
      number_of_points=len(deik[0])
      #variance[deik]=1.
      chisq=np.sum((C_out_burst[deik]-C_fit_burst[deik])**2./variance[deik])
      red_chisq=chisq/(number_of_points-4.)
      return (red_chisq)

    rr=optimize.fmin(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.],maxfun=5000)  # the initial conditions are close to the input ones, but tested to work for a broad range
    #rr=optimize.minimize(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.])
    #print (rr.x)
    n_out=rr[0]/(1e6)
    T_par_out=0.5*me*rr[1]*rr[1]/joule
    T_perp_out=0.5*me*rr[2]*rr[2]/joule
    kappa_out=rr[3]

    Density_results_array[Density_deik,sample_deik]=n_out
    Tpar_results_array[Density_deik,sample_deik]=T_par_out
    Tperp_results_array[Density_deik,sample_deik]=T_perp_out
    kappa_results_array[Density_deik,sample_deik]=kappa_out

plt.subplot(2,2,1)
plt.hist(Density_results_array[0,:])
plt.subplot(2,2,2)
plt.hist(kappa_results_array[0,:])
plt.subplot(2,2,3)
plt.hist(Tpar_results_array[0,:])
plt.subplot(2,2,4)
plt.hist(Tperp_results_array[0,:])
plt.show()

np.save('exp_counts_array.npy',exp_counts_array)
#quit()
np.save('k3_kappa_results_array.npy',kappa_results_array)      #########################################################################################
np.save('k3_Density_results_array.npy',Density_results_array)  ## Those .npy files should be stored in different folders. One for C>1 and one for C<1 ##
np.save('k3_Tpar_results_array.npy',Tpar_results_array)        ## Then the .npy files are processed by the analyze_results.py algorithm               ##
np.save('k3_Tperp_results_array.npy',Tperp_results_array)      #########################################################################################
np.save('Density_array.npy',Density_array)


print ('END')    #############
    #############






