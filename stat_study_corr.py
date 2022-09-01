#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

err_mq=20
eseed=100

path='/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_{1}seeeds_corrhplus/'.format(err_mq,eseed)

fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)

rms_dist_x_it1 =np.array(0)
ave_dist_x_it1 =np.array(0)
rms_dist_y_it1 =np.array(0)
ave_dist_y_it1 =np.array(0)
rms_dist_x_it2 =np.array(0)
ave_dist_x_it2 =np.array(0)
rms_dist_y_it2 =np.array(0)
ave_dist_y_it2 =np.array(0)

xseed = np.array(0)
iis=0

p=182 #Gev
Brho=3.3356*p

for i in range(eseed):
    
    fnameallx=path+"test_seed{0}/cx_fccee_heb_mic_all.tab".format(i+1)
    fnameally=path+"test_seed{0}/cy_fccee_heb_mic_all.tab".format(i+1)
    if os.path.exists(fnameallx):
        xseed=np.append(xseed,iis)
        head_opt=pd.read_csv(fnameallx, header=6, sep='\s+', nrows=0).columns[1:]
        optics_allx=pd.read_csv(fnameallx, skiprows=8, sep='\s+', names=head_opt)
        optics_allx=optics_allx.reset_index(drop=True)
        xcorr=np.arange(0,len(optics_allx),1)
        rms_dist_x_it1=np.append(rms_dist_x_it1,np.std(optics_allx['PX.OLD']*Brho))
        ave_dist_x_it1=np.append(ave_dist_x_it1,np.mean(optics_allx['PX.OLD']*Brho))
        rms_dist_x_it2=np.append(rms_dist_x_it2,np.std(optics_allx['PX.CORRECTION']*Brho))
        ave_dist_x_it2=np.append(ave_dist_x_it2,np.mean(optics_allx['PX.CORRECTION']*Brho))

        head_opt=pd.read_csv(fnameally, header=6, sep='\s+', nrows=0).columns[1:]
        optics_ally=pd.read_csv(fnameally, skiprows=8, sep='\s+', names=head_opt)
        optics_ally=optics_ally.reset_index(drop=True)
        ycorr=np.arange(0,len(optics_ally),1)
        rms_dist_y_it1=np.append(rms_dist_y_it1,np.std(optics_ally['PY.OLD']*Brho))
        ave_dist_y_it1=np.append(ave_dist_y_it1,np.mean(optics_ally['PY.OLD']*Brho))
        rms_dist_y_it2=np.append(rms_dist_y_it2,np.std(optics_ally['PY.CORRECTION']*Brho))
        ave_dist_y_it2=np.append(ave_dist_y_it2,np.mean(optics_ally['PY.CORRECTION']*Brho))

        
        ax3[0].plot(xcorr, optics_allx['PX.OLD']*Brho, ".")
        ax3[0].set_ylabel("Corrector strength [Tm]")
        ax3[0].set_ylim(-120e-4,120e-4)
        ax3[1].plot(ycorr, optics_ally['PY.OLD']*Brho, ".")
        ax3[1].set_xlabel("Correctors")
        ax3[1].set_ylabel("Corrector strength [Tm]")
        ax3[1].set_ylim(-120e-4,120e-4)
        fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
        
        ax5[0].plot(xcorr, optics_allx['PX.CORRECTION']*Brho, ".")
        ax5[0].set_ylabel("Corrector strength [Tm]")
        ax5[0].set_ylim(-120e-4,120e-4)
        ax5[1].plot(ycorr, optics_ally['PY.CORRECTION']*Brho, ".")
        ax5[1].set_xlabel("Correctors")
        ax5[1].set_ylabel("Corrector strength [Tm]")
        ax5[1].set_ylim(-120e-4,120e-4)
        fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
        
        iis+=1
 
#to heavy for pdf format 
fig3, plt.savefig(path+'correctors_strength_distribution_it1_{0}seeds.png'.format(eseed))
fig5, plt.savefig(path+'correctors_strength_distribution_it2_{0}seeds.png'.format(eseed))   

fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x_it1, ".")
ax2[0].plot(xseed, rms_dist_x_it2, ".")
ax2[0].set_ylabel("rms$_x$ [Tm]")
#ax2[0].set_ylim(150e-5,400e-5) #80um
ax2[0].set_ylim(0,800e-6) #20um
ax2[1].plot(xseed, rms_dist_y_it1, ".")
ax2[1].plot(xseed, rms_dist_y_it2, ".")
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [Tm]")
#ax2[1].set_ylim(150e-5,250e-5) #80um
ax2[1].set_ylim(0,800e-6) #20um
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2, plt.savefig(path+'rms_correctors_strength_{0}seeds.pdf'.format(eseed))

fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x_it1, ".")
ax4[0].plot(xseed, ave_dist_x_it2, ".")
ax4[0].set_ylabel("mean$_x$ [Tm]")
ax4[0].set_ylim(-300e-6,300e-6)
ax4[1].plot(xseed, ave_dist_y_it1, ".")
ax4[1].plot(xseed, ave_dist_y_it2, ".")
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [Tm]")
ax4[1].set_ylim(-200e-6,200e-6)
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4, plt.savefig(path+'mean_correctors_strength_{0}seeds.pdf'.format(eseed))

#for iteration n°1
fig1, ax1=plt.subplots(nrows=2, ncols=2, sharey=True)
ax1[0,0].hist(rms_dist_x_it1,100)
ax1[0,0].set_xlabel("rms$_x$ [Tm]")
ax1[0,1].hist(ave_dist_x_it1,100)
ax1[0,1].set_xlabel("mean$_x$ [Tm]")
ax1[0,0].set_ylabel("counts")
ax1[1,0].hist(rms_dist_y_it1,100)
ax1[1,0].set_xlabel("rms$_y$ [Tm]")
ax1[1,1].hist(ave_dist_y_it1,100)
ax1[1,1].set_xlabel("mean$_y$ [Tm]")
ax1[1,0].set_ylabel("counts")
fig1, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig1.tight_layout()
fig1, plt.savefig(path+'hist_correctors_strength_it1_{0}seeds.pdf'.format(eseed))

#for iteration n°2
fig8, ax8=plt.subplots(nrows=2, ncols=2, sharey=True)
ax8[0,0].hist(rms_dist_x_it2,100)
ax8[0,0].set_xlabel("rms$_x$ [Tm]")
ax8[0,1].hist(ave_dist_x_it2,100)
ax8[0,1].set_xlabel("mean$_x$ [Tm]")
ax8[0,0].set_ylabel("counts")
ax8[1,0].hist(rms_dist_y_it2,100)
ax8[1,0].set_xlabel("rms$_y$ [Tm]")
ax8[1,1].hist(ave_dist_y_it2,100)
ax8[1,1].set_xlabel("mean$_y$ [Tm]")
ax8[1,0].set_ylabel("counts")
fig8, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig8.tight_layout()
fig8, plt.savefig(path+'hist_correctors_strength_it2_{0}seeds.pdf'.format(eseed))

#plt.show()


