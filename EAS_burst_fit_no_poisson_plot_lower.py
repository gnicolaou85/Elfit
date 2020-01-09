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
######################
Density=5. * 1000000. # density in cm^-3
E0=0.7115               # (for 400 km/s E=0.455, for 500km/s E=0.712  for 800km/s E=1.82)

anisotropy=2.
Eth_par=10.
Eth_perp=anisotropy*Eth_par

uth_par=np.sqrt(2.*Eth_par*joule/me)
uth_perp=np.sqrt(2.*Eth_perp*joule/me)

kappa=3.
alpha=45.     # angle between velocity (x-axis) and B_field (on top hat plane)


print ('## INPUT PARAMETERS ##')
print ('density=',Density/(1e6))
print ('Parallel Temperature=',Eth_par,' Uth_par = ',uth_par)
print ('Perp Temperature=',Eth_perp,' Uth_perp = ',uth_perp)
print ('Speed = ',np.sqrt(2.*E0*joule/me)/1000.)
print ('kappa = ',kappa)
print ('######################')
#quit()
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

## building the INPUT distribution ##
Upar=Ux*np.cos(alpha*rad)+Uy*np.sin(alpha*rad)
Uperp= Ux*np.sin(alpha*rad)-Uy*np.cos(alpha*rad)  ### U perp on top hat plane. (Uz is also Uperp)

parenthesis=(Upar-U0_par)**2./(uth_par**2.) + (Uperp-U0_perp)**2./(uth_perp**2.)+ (Uz**2.)/(uth_perp**2.)
          
norm=np.pi**(-1.5)*(uth_par**(-1.))*(uth_perp**(-2.))*((kappa-1.5)**(-1.5))*(math.gamma(kappa+1.)/math.gamma(kappa-0.5))
f=Density*norm*((1.+parenthesis/(kappa-1.5))**(-kappa-1.))

C=f*(U**4.)*Geometric_factor*DT   ## counts
    
C_out=1.*C #np.random.poisson(C)                    # actual_measurement in 3D / here we comment the np.random.Poisson(C) if we do not wish Poisson noise
C_out_burst=C_out[:,7,:]                            # select the elevation of 0 degrees (that is the slice to me measured in burst mode for B field on the top-hat) 
#print (f.shape)
#quit()

###############################################
## plot the counts in the instrument's frame ##
###############################################
#rc('font', weight='bold')
fig=plt.figure(figsize=(6,6))
ax1 = fig.add_axes([0.15,0.25,0.6,0.6], projection='polar')
ax1.set_rticks([1, 2, 3, 4])
ax1.set_rlabel_position(20.)
ax1.tick_params(labelsize=12,labelcolor='black')
spec=ax1.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(C_out_burst),cmap=cmap)
cax1=fig.add_axes([0.15,0.13,0.6,0.04])
cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(Counts)',cax=cax1)#,pad=0.08,fraction=0.054)
ax1.set_xlim([0.*rad,360.*rad])
ax1.grid(linestyle='--',alpha=0.5)
plt.show()
###############################################
###############################################

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

#####################################################
## plot the distribution in the instrument's frame ##
#####################################################
#fig=plt.figure(figsize=(6,6))
#ax1 = fig.add_axes([0.15,0.25,0.6,0.6], projection='polar')
#ax1.set_rticks([1, 2, 3, 4])
#ax1.set_rlabel_position(20.)
#ax1.tick_params(labelsize=12,labelcolor='black')
#spec=ax1.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(f_out_burst),cmap=cmap)
#cax1=fig.add_axes([0.15,0.13,0.6,0.04])
#cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(f_out)',cax=cax1)#,pad=0.08,fraction=0.054)
#ax1.set_xlim([0.*rad,360.*rad])
#ax1.grid(linestyle='--',alpha=0.5)
#plt.show()
#####################################################
#####################################################



##################
## Fitting C>=1 ##
##################
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
  #plt.show()
  variance=1.*C_out_burst
  deik=np.where(C_out_burst>=1.)
  number_of_points=len(deik[0]) #len(C_out_burst)*len(C_out_burst[0])
  #variance[deik]=1.
  chisq=np.sum((C_out_burst[deik]-C_fit_burst[deik])**2./variance[deik])
  red_chisq=chisq/(number_of_points-4.)
  print (number_of_points,red_chisq)
  return (red_chisq)

rr=optimize.fmin(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.],maxfun=1000) # initial guessings are chosen close to the input values, but tested for a wide range of values
#rr=optimize.minimize(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.])

#print (rr)
n_out=rr[0]/(1e6)
T_par_out=0.5*me*rr[1]*rr[1]/joule
T_perp_out=0.5*me*rr[2]*rr[2]/joule
kappa_out=rr[3]

print (n_out,T_par_out,T_perp_out,kappa_out)


##########################
## plot fitted function ##
##########################


parf=(Upar-U0_par)**2./(rr[1]**2.) + (Uperp-U0_perp)**2./(rr[2]**2.)+ (Uz**2.)/(rr[2]**2.)
          
normf=np.pi**(-1.5)*(rr[1]**(-1.))*(rr[2]**(-2.))*((rr[3]-1.5)**(-1.5))*(math.gamma(rr[3]+1.)/math.gamma(rr[3]-0.5))
f_fit=rr[0]*normf*((1.+parf/(rr[3]-1.5))**(-rr[3]-1.))

C_fit=f_fit*(U**4.)*Geometric_factor*DT

C_fit_burst=C_fit[:,7,:]

##########################################################################
## plot the counts in the instrument's frame and the fitted model aside ##
##########################################################################
#rc('font', weight='bold')
fig=plt.figure(figsize=(13,6))
ax1 = fig.add_axes([0.15,0.25,0.3,0.6], projection='polar')
ax1.set_rticks([1, 2, 3, 4])
ax1.set_rlabel_position(20.)
ax1.tick_params(labelsize=12,labelcolor='black')
spec=ax1.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(C_out_burst),cmap=cmap,vmin=0,vmax=np.max(np.log10(C_out_burst)))
cax1=fig.add_axes([0.35,0.13,0.3,0.04])
cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(Counts)',cax=cax1)#,pad=0.08,fraction=0.054)
ax1.set_xlim([0.*rad,360.*rad])
ax1.grid(linestyle='--',alpha=0.5)

ax2 = fig.add_axes([0.55,0.25,0.3,0.6], projection='polar')
ax2.set_rticks([1, 2, 3, 4])
ax2.set_rlabel_position(20.)
ax2.tick_params(labelsize=12,labelcolor='black')
spec=ax2.pcolor(azimuth_plot*rad,np.log10(energy_plot),np.log10(C_fit_burst),cmap=cmap,vmin=0,vmax=np.max(np.log10(C_out_burst)))
#cax2=fig.add_axes([0.55,0.13,0.3,0.04])
#cbar=plt.colorbar(spec,extend='min',orientation='horizontal',label=r'log$_{10}$(Counts)',cax=cax2)#,pad=0.08,fraction=0.054)
ax2.set_xlim([0.*rad,360.*rad])
ax2.grid(linestyle='--',alpha=0.5)
plt.show()
#quit()
#########################################################################
#########################################################################

line_plot=C_out_burst[:,28]
fit_plot=C_fit_burst[:,28]

i=np.where(line_plot>=1.)
j=np.where(line_plot<1.)

fig=plt.figure(figsize=(8,4))
ax1 = fig.add_axes([0.1,0.14,0.42,0.75])
ax1.errorbar(energy_table_middle[i],line_plot[i],yerr=np.sqrt(line_plot[i]),elinewidth=1,capsize=2,fmt='o',color='black',label='measured',zorder=-1)
ax1.errorbar(energy_table_middle[j],line_plot[j],yerr=0,color='red',fmt='o',zorder=-1)
ax1.plot(energy_table_middle,fit_plot,color='blue',alpha=1.,linewidth=2,label='fitted')
ax1.plot(energy_table_middle,C[:,7,28]-1.,color='orange',linewidth=2,alpha=1,label=r'$C_{\mathrm{exp}}$-1')
ax1.set_xscale('log')
ax1.tick_params(labelsize=12)
ax1.set_ylabel('Counts',fontsize=14)
ax1.set_xlabel('Energy (eV)',fontsize=14)
ax1.set_ylim([0,15])#1.2*np.max(C_out_burst)])
ax1.tick_params(labelsize=12)
ax1.grid(linestyle='--',alpha=0.5)
ax1.text(200,10,r'$n$ = '+str(round(n_out,1))+r' cm$^{-3}$')
ax1.text(200,9,r'$T_{\parallel}$ = '+str(round(T_par_out,1))+r' eV')
ax1.text(200,8,r'$T_{\perp}$ = '+str(round(T_perp_out,1))+r' eV')
ax1.text(200,7,r'$\kappa_{\mathrm{out}}$ = '+str(np.round(kappa_out,1)))
ax1.text(200,14,r'fit $C_{i} \geq$ 1',bbox=dict(facecolor='green', alpha=0.5))
ax1.legend(loc=2)

#plt.show()
#############
#############

#####################
## Fitting every C ##
#####################
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
  #plt.show()
  variance=1.*C_out_burst
  deik=np.where(C_out_burst<1.)
  number_of_points=len(C_out_burst)*len(C_out_burst[0])
  variance[deik]=1.
  chisq=np.sum((C_out_burst-C_fit_burst)**2./variance)
  red_chisq=chisq/(number_of_points-4.)
  print (number_of_points,red_chisq)
  return (red_chisq)

rr=optimize.fmin(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.],maxfun=1000)    # initial guessings are chosen close to the input values, but tested for a wide range of values
#rr=optimize.minimize(fopt, [Density,np.sqrt(2.*10.*joule/me),np.sqrt(2.*10.*joule/me),5.])

#print (rr)
n_out=rr[0]/(1e6)
T_par_out=0.5*me*rr[1]*rr[1]/joule
T_perp_out=0.5*me*rr[2]*rr[2]/joule
kappa_out=rr[3]

print (n_out,T_par_out,T_perp_out,kappa_out)


##################################
## plot fitted & input function ##
##################################


parf=(Upar-U0_par)**2./(rr[1]**2.) + (Uperp-U0_perp)**2./(rr[2]**2.)+ (Uz**2.)/(rr[2]**2.)
          
normf=np.pi**(-1.5)*(rr[1]**(-1.))*(rr[2]**(-2.))*((rr[3]-1.5)**(-1.5))*(math.gamma(rr[3]+1.)/math.gamma(rr[3]-0.5))
f_fit=rr[0]*normf*((1.+parf/(rr[3]-1.5))**(-rr[3]-1.))

C_fit=f_fit*(U**4.)*Geometric_factor*DT

C_fit_burst=C_fit[:,7,:]


line_plot=C_out_burst[:,28]
fit_plot=C_fit_burst[:,28]

i=np.where(line_plot>=1.)
j=np.where(line_plot<1.)

ax2 = fig.add_axes([0.54,0.14,0.42,0.75])
ax2.errorbar(energy_table_middle[i],line_plot[i],yerr=np.sqrt(line_plot[i]),elinewidth=1,capsize=2,fmt=' o',color='black',label='measured',zorder=-1)
ax2.errorbar(energy_table_middle[j],line_plot[j],yerr=1,elinewidth=1,capsize=2,color='black',fmt=' o',zorder=-1)
ax2.plot(energy_table_middle,fit_plot,color='blue',linewidth=2,label='fitted')
ax2.plot(energy_table_middle,C[:,7,28]-1.,color='orange',linewidth=2,alpha=1,label=r'$C_{\mathrm{exp}}$-1')
ax2.set_xscale('log')
ax2.set_ylim([0,15])#1.2*np.max(C_out_burst)])
ax2.tick_params(labelsize=12,labelleft=False)
ax2.set_xlabel('Energy (eV)',fontsize=14)
ax2.grid(linestyle='--',alpha=0.5)
ax2.text(200,10,r'$n$ = '+str(round(n_out,1))+r' cm$^{-3}$')
ax2.text(200,9,r'$T_{\parallel}$ = '+str(round(T_par_out,1))+r' eV')
ax2.text(200,8,r'$T_{\perp}$ = '+str(round(T_perp_out,1))+r' eV')
ax2.text(200,7,r'$\kappa_{\mathrm{out}}$ = '+str(np.round(kappa_out,1)))
ax2.text(200,14,r'fit all',bbox=dict(facecolor='green', alpha=0.5))
plt.title(r'Input $n\,=\,5\,$cm$^{-3}$,  $T_{\parallel}\,=\,10\,$eV,  $T_{\perp}\,=\,20\,$eV,  $\kappa\,=\,3$',y=1.015,x=-0.1)
ax2.legend(loc=2)
plt.show()
#############
#############


