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

rms_dist_x =np.array(0)
ave_dist_x =np.array(0)
rms_dist_y =np.array(0)
ave_dist_y =np.array(0)
rms_dist_x2 =np.array(0)
ave_dist_x2 =np.array(0)
rms_dist_y2 =np.array(0)
ave_dist_y2 =np.array(0)

xseed = np.array(0)
xseed2 = np.array(0)
iis=0
jjs=0

for i in range(eseed):
    
    #for iteration n°1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall1):
        xseed=np.append(xseed,iis)
        head_opt=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt)
        optics_all1=optics_all1.reset_index(drop=True)
        rms_dist_x=np.append(rms_dist_x,np.std(optics_all1['X']))
        ave_dist_x=np.append(ave_dist_x,np.mean(optics_all1['X']))
        rms_dist_y=np.append(rms_dist_y,np.std(optics_all1['Y']))
        ave_dist_y=np.append(ave_dist_y,np.mean(optics_all1['Y']))

        ax3[0].plot(optics_all1['S']/1000., optics_all1['X'], ".")
        ax3[0].set_ylabel("x [m]")
        ax3[0].set_ylim(-300e-5,300e-5)
        ax3[1].plot(optics_all1['S']/1000., optics_all1['Y'], ".")
        ax3[1].set_xlabel("longitundinal position [km]")
        ax3[1].set_ylabel("y [m]")
        ax3[1].set_ylim(-500e-6,500e-6)
        fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

        iis+=1
    else:
        continue
    
    #for iteration n°2
    fnameall2=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall2):
        xseed2=np.append(xseed2,jjs)
        head_opt=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all2=pd.read_csv(fnameall2, skiprows=52, sep='\s+', names=head_opt)
        optics_all2=optics_all2.reset_index(drop=True)
        rms_dist_x2=np.append(rms_dist_x2,np.std(optics_all2['X']))
        ave_dist_x2=np.append(ave_dist_x2,np.mean(optics_all2['X']))
        rms_dist_y2=np.append(rms_dist_y2,np.std(optics_all2['Y']))
        ave_dist_y2=np.append(ave_dist_y2,np.mean(optics_all2['Y']))

        ax5[0].plot(optics_all2['S']/1000., optics_all2['X'], ".")
        ax5[0].set_ylabel("x [m]")
        ax5[0].set_ylim(-300e-5,300e-5)
        ax5[1].plot(optics_all2['S']/1000., optics_all2['Y'], ".")
        ax5[1].set_xlabel("longitundinal position [km]")
        ax5[1].set_ylabel("y [m]")
        ax5[1].set_ylim(-150e-6,150e-6)
        fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
 
        jjs+=1
    else:
        continue

#to heavy for pdf format    
fig3, plt.savefig(path+'orbit_distribution_it1_{0}seeds.png'.format(eseed))
fig5, plt.savefig(path+'orbit_distribution_it2_{0}seeds.png'.format(eseed)) 

fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x, ".")
ax2[0].plot(xseed2, rms_dist_x2, ".")
ax2[0].set_ylabel("rms$_x$ [m]")
ax2[0].set_ylim(-50e-6,200e-5)
ax2[1].plot(xseed, rms_dist_y, ".")
ax2[1].plot(xseed2, rms_dist_y2, ".")
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m]")
ax2[1].set_ylim(0,500e-6)
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2, plt.savefig(path+'rms_orbit_{0}seeds.pdf'.format(eseed))

fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x, ".")
ax4[0].plot(xseed2, ave_dist_x2, ".")
ax4[0].set_ylabel("mean$_x$ [m]")
ax4[0].set_ylim(-50e-6,50e-6)
ax4[1].plot(xseed, ave_dist_y, ".")
ax4[1].plot(xseed2, ave_dist_y2, ".")
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m]")
ax4[1].set_xlim(0,93)
ax4[1].set_ylim(-50e-6,50e-6)
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4, plt.savefig(path+'mean_orbit_{0}seeds.pdf'.format(eseed))

#for iteration n°1
fig1, ax1=plt.subplots(nrows=2, ncols=2, sharey=True)
ax1[0,0].hist(rms_dist_x,100)
ax1[0,0].set_xlabel("rms$_x$ [m]")
ax1[0,1].hist(ave_dist_x,100)
ax1[0,1].set_xlabel("mean$_x$ [m]")
ax1[0,0].set_ylabel("counts")
ax1[1,0].hist(rms_dist_y,100)
ax1[1,0].set_xlabel("rms$_y$ [m]")
ax1[1,1].hist(ave_dist_y,100)
ax1[1,1].set_xlabel("mean$_y$ [m]")
ax1[1,0].set_ylabel("counts")
fig1, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig1.tight_layout()
fig1, plt.savefig(path+'hist_orbit_it1_{0}seeds.pdf'.format(eseed))

#for iteration n°2
fig8, ax8=plt.subplots(nrows=2, ncols=2, sharey=True)
ax8[0,0].hist(rms_dist_x2,100)
ax8[0,0].set_xlabel("rms$_x$ [m]")
ax8[0,1].hist(ave_dist_x2,100)
ax8[0,1].set_xlabel("mean$_x$ [m]")
ax8[0,0].set_ylabel("counts")
ax8[1,0].hist(rms_dist_y2,100)
ax8[1,0].set_xlabel("rms$_y$ [m]")
ax8[1,1].hist(ave_dist_y2,100)
ax8[1,1].set_xlabel("mean$_y$ [m]")
ax8[1,0].set_ylabel("counts")
fig8, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig8.tight_layout()
fig8, plt.savefig(path+'hist_orbit_it2_{0}seeds.pdf'.format(eseed))

#plt.show()


