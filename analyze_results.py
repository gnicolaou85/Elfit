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

k=np.load('counts_gt_one/k3_kappa_results_array.npy')                     #########################################
n=np.load('counts_gt_one/k3_Density_results_array.npy')                   ##the name of the folders to read from ##
Tpar=np.load('counts_gt_one/k3_Tpar_results_array.npy')                   #########################################
Tperp=np.load('counts_gt_one/k3_Tperp_results_array.npy')
n_in=np.load('counts_gt_one/Density_array.npy')
exp_counts_array=np.load('counts_gt_one/exp_counts_array.npy')

n_mean=np.zeros(len(n_in))
n_std=np.zeros(len(n_in))
n_median=np.zeros(len(n_in))
n_lowerq=np.zeros(len(n_in))
n_higherq=np.zeros(len(n_in))


k_mean=np.zeros(len(n_in))
k_std=np.zeros(len(n_in))
k_median=np.zeros(len(n_in))
k_higherq=np.zeros(len(n_in))
k_lowerq=np.zeros(len(n_in))


Tpar_mean=np.zeros(len(n_in))
Tpar_std=np.zeros(len(n_in))
Tpar_median=np.zeros(len(n_in))
Tpar_higherq=np.zeros(len(n_in))
Tpar_lowerq=np.zeros(len(n_in))

Tperp_mean=np.zeros(len(n_in))
Tperp_std=np.zeros(len(n_in))
Tperp_median=np.zeros(len(n_in))
Tperp_higherq=np.zeros(len(n_in))
Tperp_lowerq=np.zeros(len(n_in))


for i in range(len(n_in)):

  n_mean[i]=np.mean(n[i,:])
  n_std[i]=np.std(n[i,:])
  n_median[i]=np.percentile(n[i,:], 50)
  n_lowerq[i]=np.percentile(n[i,:], 25)
  n_higherq[i]=np.percentile(n[i,:], 75)

  k_mean[i]=np.mean(k[i,:])
  k_std[i]=np.std(k[i,:])
  k_median[i]=np.percentile(k[i,:], 50)
  k_lowerq[i]=np.percentile(k[i,:], 25)
  k_higherq[i]=np.percentile(k[i,:], 75)

  Tpar_mean[i]=np.mean(Tpar[i,:])
  Tpar_std[i]=np.std(Tpar[i,:])
  Tpar_median[i]=np.percentile(Tpar[i,:], 50)
  Tpar_lowerq[i]=np.percentile(Tpar[i,:], 25)
  Tpar_higherq[i]=np.percentile(Tpar[i,:], 75)

  Tperp_mean[i]=np.mean(Tperp[i,:])
  Tperp_std[i]=np.std(Tperp[i,:])
  Tperp_median[i]=np.percentile(Tperp[i,:], 50)
  Tperp_lowerq[i]=np.percentile(Tperp[i,:], 25)
  Tperp_higherq[i]=np.percentile(Tperp[i,:], 75)


col1='red'


fig=plt.figure(figsize=(7,7))
ax1 = fig.add_axes([0.1,0.70,0.8,0.18])
ax1.plot(n_in*(10.**(-6.)),n_in*(10.**(-6.))*0.+1.,'--',color='black',linewidth=2)
ax1.plot(n_in*(10.**(-6.)),n_mean/(n_in*(10.**(-6.))),'o',color=col1,label=r'fit $C \geq$ 1')
ax1.fill_between(n_in*(10.**(-6.)),(n_mean-n_std)/(n_in*(10.**(-6.))),(n_mean+n_std)/(n_in*(10.**(-6.))),facecolor=col1,alpha=0.3)
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.set_ylim([0.5,1.5])
ax1.grid(linestyle='--',alpha=0.5)

ax11 = ax1.twiny()
ax11.plot(exp_counts_array,n_mean/(n_in*(10.**(-6.))),'o',color=col1) # Create a dummy plot
ax11.set_xscale('log')
#ax11.set_yscale('log')
ax11.set_ylim([0.5,1.5])

ax11.set_xlabel('max Counts',color='grey',fontsize=12)
ax11.tick_params(labelsize=12,colors='grey')

ax2 = fig.add_axes([0.1,0.50,0.8,0.18])
ax2.plot(n_in*(10.**(-6.)),np.zeros(len(n_in))+3.,'--',color='black',linewidth=2)
ax2.plot(n_in*(10.**(-6.)),k_mean,'o',color=col1)
ax2.fill_between(n_in*(10.**(-6.)),k_mean-k_std,k_mean+k_std,facecolor=col1,alpha=0.3)
ax2.set_xscale('log')
ax2.set_ylim([1.5,4])
ax2.grid(linestyle='--',alpha=0.5)

ax3 = fig.add_axes([0.1,0.30,0.8,0.18])
ax3.plot(n_in*(10.**(-6.)),np.zeros(len(n_in))+10.,'--',color='black',linewidth=2,label='input')
ax3.plot(n_in*(10.**(-6.)),Tpar_mean,'o',color=col1,label=r'exclude $C_{i}$ = 0')
ax3.fill_between(n_in*(10.**(-6.)),Tpar_mean-Tpar_std,Tpar_mean+Tpar_std,facecolor=col1,alpha=0.3)
ax3.set_xscale('log')
ax3.set_ylim([5,25])
ax3.grid(linestyle='--',alpha=0.5)

ax4 = fig.add_axes([0.1,0.1,0.8,0.18])
ax4.plot(n_in*(10.**(-6.)),np.zeros(len(n_in))+20.,'--',color='black',linewidth=2)
ax4.plot(n_in*(10.**(-6.)),Tperp_mean,'o',color=col1)
ax4.fill_between(n_in*(10.**(-6.)),Tperp_mean-Tperp_std,Tperp_mean+Tperp_std,facecolor=col1,alpha=0.3)
ax4.set_xscale('log')
ax4.set_ylim([15,35])
ax4.grid(linestyle='--',alpha=0.5)


#print (n_mean*(1e6)/n_in)
#quit()
#plt.show()
##################
### C<1 fitting ##
##################

k=np.load('counts_lt_one/k3_kappa_results_array.npy')
n=np.load('counts_lt_one/k3_Density_results_array.npy')
Tpar=np.load('counts_lt_one/k3_Tpar_results_array.npy')
Tperp=np.load('counts_lt_one/k3_Tperp_results_array.npy')
n_in=np.load('counts_lt_one/Density_array.npy')

n_mean=np.zeros(len(n_in))
n_std=np.zeros(len(n_in))
n_median=np.zeros(len(n_in))
n_lowerq=np.zeros(len(n_in))
n_higherq=np.zeros(len(n_in))


k_mean=np.zeros(len(n_in))
k_std=np.zeros(len(n_in))
k_median=np.zeros(len(n_in))
k_higherq=np.zeros(len(n_in))
k_lowerq=np.zeros(len(n_in))


Tpar_mean=np.zeros(len(n_in))
Tpar_std=np.zeros(len(n_in))
Tpar_median=np.zeros(len(n_in))
Tpar_higherq=np.zeros(len(n_in))
Tpar_lowerq=np.zeros(len(n_in))

Tperp_mean=np.zeros(len(n_in))
Tperp_std=np.zeros(len(n_in))
Tperp_median=np.zeros(len(n_in))
Tperp_higherq=np.zeros(len(n_in))
Tperp_lowerq=np.zeros(len(n_in))


for i in range(len(n_in)):

  n_mean[i]=np.mean(n[i,:])
  n_std[i]=np.std(n[i,:])
  n_median[i]=np.percentile(n[i,:], 50)
  n_lowerq[i]=np.percentile(n[i,:], 25)
  n_higherq[i]=np.percentile(n[i,:], 75)

  k_mean[i]=np.mean(k[i,:])
  k_std[i]=np.std(k[i,:])
  k_median[i]=np.percentile(k[i,:], 50)
  k_lowerq[i]=np.percentile(k[i,:], 25)
  k_higherq[i]=np.percentile(k[i,:], 75)

  Tpar_mean[i]=np.mean(Tpar[i,:])
  Tpar_std[i]=np.std(Tpar[i,:])
  Tpar_median[i]=np.percentile(Tpar[i,:], 50)
  Tpar_lowerq[i]=np.percentile(Tpar[i,:], 25)
  Tpar_higherq[i]=np.percentile(Tpar[i,:], 75)

  Tperp_mean[i]=np.mean(Tperp[i,:])
  Tperp_std[i]=np.std(Tperp[i,:])
  Tperp_median[i]=np.percentile(Tperp[i,:], 50)
  Tperp_lowerq[i]=np.percentile(Tperp[i,:], 25)
  Tperp_higherq[i]=np.percentile(Tperp[i,:], 75)


col1='blue'
fontsize=12
labelsize=10

#fig=plt.figure(figsize=(7,7))
#ax1 = fig.add_axes([0.1,0.76,0.8,0.2])
ax1.plot(n_in*(10.**(-6.)),n_mean/(n_in*(10.**(-6.))),'o',color=col1)
ax1.fill_between(n_in*(10.**(-6.)),(n_mean-n_std)/(n_in*(10.**(-6.))),(n_mean+n_std)/(n_in*(10.**(-6.))),facecolor=col1,alpha=0.3)
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.tick_params(labelsize=labelsize,labelbottom=False)
ax1.set_ylabel('$n_{\mathrm{out}}$ / $n$',fontsize=fontsize)
#ax1.grid(linestyle='--',alpha=0.5)
ax1.set_ylim([0.5,1.5])
ax11.plot(exp_counts_array,n_mean/(n_in*(10.**(-6.))),'o',color=col1) # Create a dummy plot
ax11.set_ylim([0.5,1.5])


#ax2 = fig.add_axes([0.1,0.54,0.8,0.2])
ax2.plot(n_in*(10.**(-6.)),k_mean,'o',color=col1)
ax2.fill_between(n_in*(10.**(-6.)),k_mean-k_std,k_mean+k_std,facecolor=col1,alpha=0.3)
ax2.set_xscale('log')
ax2.tick_params(labelsize=labelsize,labelbottom=False)
ax2.set_ylabel('$\kappa_{\mathrm{out}}$',fontsize=fontsize)
#ax2.grid(linestyle='--',alpha=0.5)

#ax3 = fig.add_axes([0.1,0.32,0.8,0.2])
ax3.plot(n_in*(10.**(-6.)),Tpar_mean,'o',color=col1,label=r'include $C_{i}$ = 0')
ax3.fill_between(n_in*(10.**(-6.)),Tpar_mean-Tpar_std,Tpar_mean+Tpar_std,facecolor=col1,alpha=0.3)
ax3.set_xscale('log')
ax3.tick_params(labelsize=labelsize,labelbottom=False)
ax3.set_ylabel('$T_{\parallel,\mathrm{out}}$ (eV)',fontsize=fontsize)
ax3.legend(bbox_to_anchor=(1., 1.3),fontsize=fontsize,shadow=True,fancybox=True)

#ax3.grid(linestyle='--',alpha=0.5)

#ax4 = fig.add_axes([0.1,0.1,0.8,0.2])
ax4.plot(n_in*(10.**(-6.)),Tperp_mean,'o',color=col1)
ax4.fill_between(n_in*(10.**(-6.)),Tperp_mean-Tperp_std,Tperp_mean+Tperp_std,facecolor=col1,alpha=0.3)
ax4.set_xscale('log')
ax4.tick_params(labelsize=labelsize)
ax4.set_ylabel('$T_{\perp,\mathrm{out}}$ (eV)',fontsize=fontsize)
ax4.set_xlabel(r'Density (cm$^{-3}$)',fontsize=fontsize)
#ax4.grid(linestyle='--',alpha=0.5)
plt.savefig('derived_vs_n.pdf')
plt.show()

#print (n_mean*(1e6)/n_in)
