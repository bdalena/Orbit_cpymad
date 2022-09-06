#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the orbit after correction.
'''

import matplotlib.pyplot as plt
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os

err_mq=90 #quads offset
eseed=100 #nb of seeds

path='/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_100seeeds_corrhplus/'.format(err_mq)

fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

rms_dist_x =np.empty(0)
ave_dist_x =np.empty(0)
rms_dist_y =np.empty(0)
ave_dist_y =np.empty(0)
rms_dist_x2 =np.empty(0)
ave_dist_x2 =np.empty(0)
rms_dist_y2 =np.empty(0)
ave_dist_y2 =np.empty(0)

xseed = np.empty(0)
xseed2 = np.empty(0)
iis=0
jjs=0

orbit_x_all_seed=np.empty(0)
orbit_y_all_seed=np.empty(0)

for i in range(eseed):

    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(i+1)
    
    mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)
    orbit_x_all_seed=np.append(orbit_x_all_seed,orbit_x)
    orbit_y_all_seed=np.append(orbit_y_all_seed,orbit_y)

    #for iteration n°1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall1):
        xseed=np.append(xseed,iis)
        head_opt1=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt1)
        optics_all1=optics_all1.reset_index(drop=True)
        rms_dist_x=np.append(rms_dist_x,np.std(optics_all1['X']))
        ave_dist_x=np.append(ave_dist_x,np.mean(optics_all1['X']))
        rms_dist_y=np.append(rms_dist_y,np.std(optics_all1['Y']))
        ave_dist_y=np.append(ave_dist_y,np.mean(optics_all1['Y']))

        ax3[0].plot(optics_all1['S']/1000., optics_all1['X'], ".")
        ax3[0].set_ylabel("x [m]")
        #ax3[0].set_ylim(-300e-6,300e-6) #for 80um
        ax3[0].set_ylim(-500e-6,500e-6) #for 100um
        ax3[1].plot(optics_all1['S']/1000., optics_all1['Y'], ".")
        ax3[1].set_xlabel("longitundinal position [km]")
        ax3[1].set_ylabel("y [m]")
        #ax3[1].set_ylim(-300e-6,300e-6) #for 80um
        ax3[1].set_ylim(-500e-6,500e-6) #for 100um
        
        iis+=1
    else:
        continue
    
    #for iteration n°2
    fnameall2=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall2):
        xseed2=np.append(xseed2,jjs)
        head_opt2=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all2=pd.read_csv(fnameall2, skiprows=52, sep='\s+', names=head_opt2)
        optics_all2=optics_all2.reset_index(drop=True)
        rms_dist_x2=np.append(rms_dist_x2,np.std(optics_all2['X']))
        ave_dist_x2=np.append(ave_dist_x2,np.mean(optics_all2['X']))
        rms_dist_y2=np.append(rms_dist_y2,np.std(optics_all2['Y']))
        ave_dist_y2=np.append(ave_dist_y2,np.mean(optics_all2['Y']))

        ax5[0].plot(optics_all2['S']/1000., optics_all2['X'], ".")
        ax5[0].set_ylabel("x [m]")
        ax5[0].set_ylim(-150e-6,150e-6) #not for 20 & 40um
        ax5[1].plot(optics_all2['S']/1000., optics_all2['Y'], ".")
        ax5[1].set_xlabel("longitundinal position [km]")
        ax5[1].set_ylabel("y [m]")
        ax5[1].set_ylim(-150e-6,150e-6) #not for 20 & 40um
        
        jjs+=1
    else:
        continue


max_orbit_x=np.max(orbit_x_all_seed)
max_orbit_y=np.max(orbit_y_all_seed)

#to heavy for pdf format
#fig3.savefig(path+'orbit_distribution_it1_{0}seeds.png'.format(eseed))
#fig5.savefig(path+'orbit_distribution_it2_{0}seeds.png'.format(eseed)) 

#plot of the rms
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x, ".")
ax2[0].plot(xseed2, rms_dist_x2, ".")
#ax2[0].axhline(y=min_orbit_x, color='g', linestyle='--')
ax2[0].axhline(y=max_orbit_x, color='r', linestyle='--')
ax2[0].set_ylabel("rms$_x$ [m]")
#ax2[0].set_ylim(0,200e-5) #for 80um
#ax2[0].set_ylim(0,6e-6) #for 20um
#ax2[0].set_ylim(0,1.9e-5) #for 40um
#ax2[0].set_ylim(0,5e-5) #for 60um
ax2[0].set_ylim(0,max_orbit_x+15e-5)
ax2[1].plot(xseed, rms_dist_y, ".")
ax2[1].plot(xseed2, rms_dist_y2, ".")
#ax2[1].axhline(y=min_orbit_y, color='g', linestyle='--')
ax2[1].axhline(y=max_orbit_y, color='r', linestyle='--')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m]")
#ax2[1].set_ylim(0,500e-6) #for 80um
#ax2[1].set_ylim(4e-6,6e-6) #for 20um
#ax2[1].set_ylim(0.9e-5,1.9e-5) #for 40um
#ax2[1].set_ylim(1e-5,5e-5) #for 60um
ax2[1].set_ylim(0,max_orbit_y+15e-5)
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_orbit_{0}seeds.pdf'.format(eseed))
sys.exit()
#plot of the mean
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x, ".")
ax4[0].plot(xseed2, ave_dist_x2, ".")
ax4[0].set_ylabel("mean$_x$ [m]")
ax4[0].set_ylim(-50e-6,50e-6) #for 80um
#ax4[0].set_ylim(-3e-7,3e-7) #for 20um
#ax4[0].set_ylim(-6e-7,6e-7) #for 40um
#ax4[0].set_ylim(-15e-7,15e-7) #for 60um
ax4[1].plot(xseed, ave_dist_y, ".")
ax4[1].plot(xseed2, ave_dist_y2, ".")
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m]")
ax4[1].set_ylim(-50e-6,50e-6) #for 80um
#ax4[1].set_ylim(-3e-7,3e-7) #for 20um
#ax4[1].set_ylim(-6e-7,6e-7) #for 40um
#ax4[1].set_ylim(-15e-7,15e-7) #for 60um
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_orbit_{0}seeds.pdf'.format(eseed))

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
fig1.savefig(path+'hist_orbit_it1_{0}seeds.pdf'.format(eseed))

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
fig8.savefig(path+'hist_orbit_it2_{0}seeds.pdf'.format(eseed))

#plt.show()


