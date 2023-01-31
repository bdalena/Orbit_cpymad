#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the orbit after correction and tune match.
'''

import matplotlib.pyplot as plt
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os

err_mq=100 #quads offset
eseed=100 #nb of seeds

#path definition
#path='/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_100seeeds_corrhplus/'.format(err_mq)
#path='/home/td271008/work/cpymadx/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_2it_100seeeds_corrhplus/'.format(err_mq)
path='./'

fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig9, ax9=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig9, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig10, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

#CORRECT: it 1
rms_dist_x =np.empty(0)
ave_dist_x =np.empty(0)
rms_dist_y =np.empty(0)
ave_dist_y =np.empty(0)

#CORRECT: it 2
rms_dist_x2 =np.empty(0)
ave_dist_x2 =np.empty(0)
rms_dist_y2 =np.empty(0)
ave_dist_y2 =np.empty(0)

#tune match: first it
rms_dist_x3 =np.empty(0)
ave_dist_x3 =np.empty(0)
rms_dist_y3 =np.empty(0)
ave_dist_y3 =np.empty(0)

#tune match: last it
rms_dist_x4 =np.empty(0)
ave_dist_x4 =np.empty(0)
rms_dist_y4 =np.empty(0)
ave_dist_y4 =np.empty(0)

xseed = np.empty(0)
xseed2 = np.empty(0)
xseed3 = np.empty(0)
xseed4 = np.empty(0)

iis=0
jjs=0
kks=0
pps=0

orbit_x_all_seed=np.empty(0)
orbit_y_all_seed=np.empty(0)  

for i in range(eseed):
    iseed= i+1
    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(iseed)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(iseed)
    
    mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)
    orbit_x_all_seed=np.append(orbit_x_all_seed,orbit_x)
    orbit_y_all_seed=np.append(orbit_y_all_seed,orbit_y)

    #CORRECT: it 1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(iseed)
    if os.path.exists(fnameall1):
        xseed=np.append(xseed,iis)
        head_opt1=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt1)
        optics_all1=optics_all1.reset_index(drop=True)
        rms_dist_x=np.append(rms_dist_x,np.std(optics_all1['X']))
        ave_dist_x=np.append(ave_dist_x,np.mean(optics_all1['X']))
        rms_dist_y=np.append(rms_dist_y,np.std(optics_all1['Y']))
        ave_dist_y=np.append(ave_dist_y,np.mean(optics_all1['Y']))
        print("correct sex on it 1= {0}".format(np.std(optics_all1['X'])))
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
    
    #CORRECT: it 2
    fnameall2=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(iseed)
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

    #tune match: first it
    fnameall3=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs".format(iseed)
    if os.path.exists(fnameall3):
        xseed3=np.append(xseed3,kks)
        head_opt3=pd.read_csv(fnameall3, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all3=pd.read_csv(fnameall3, skiprows=52, sep='\s+', names=head_opt3)
        optics_all3=optics_all3.reset_index(drop=True)
        rms_dist_x3=np.append(rms_dist_x3,np.std(optics_all3['X']))
        ave_dist_x3=np.append(ave_dist_x3,np.mean(optics_all3['X']))
        rms_dist_y3=np.append(rms_dist_y3,np.std(optics_all3['Y']))
        ave_dist_y3=np.append(ave_dist_y3,np.mean(optics_all3['Y']))

        ax9[0].plot(optics_all3['S']/1000., optics_all3['X'], ".")
        ax9[0].set_ylabel("x [m]")
        ax9[0].set_ylim(-4e-4,4e-4) #not for 20 & 40um
        ax9[1].plot(optics_all3['S']/1000., optics_all3['Y'], ".")
        ax9[1].set_xlabel("longitundinal position [km]")
        ax9[1].set_ylabel("y [m]")
        ax9[1].set_ylim(-2e-4,2e-4) #not for 20 & 40um
        
        kks+=1
    else:
        continue

    cnt=1
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(iseed,cnt)
    
    while os.path.exists(fnameall4):
        fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(iseed,cnt+1)
        cnt+=1
    
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(iseed,cnt-1)
    
    #tune match: last it
    if os.path.exists(fnameall4):
        xseed4=np.append(xseed4,pps)
        head_opt4=pd.read_csv(fnameall4, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all4=pd.read_csv(fnameall4, skiprows=52, sep='\s+', names=head_opt4)
        optics_all4=optics_all4.reset_index(drop=True)
        rms_dist_x4=np.append(rms_dist_x4,np.std(optics_all4['X']))
        ave_dist_x4=np.append(ave_dist_x4,np.mean(optics_all4['X']))
        rms_dist_y4=np.append(rms_dist_y4,np.std(optics_all4['Y']))
        ave_dist_y4=np.append(ave_dist_y4,np.mean(optics_all4['Y']))

        ax10[0].plot(optics_all4['S']/1000., optics_all4['X'], ".")
        ax10[0].set_ylabel("x [m]")
        ax10[0].set_ylim(-2e-4,2e-4) #not for 20 & 40um
        ax10[1].plot(optics_all4['S']/1000., optics_all4['Y'], ".")
        ax10[1].set_xlabel("longitundinal position [km]")
        ax10[1].set_ylabel("y [m]")
        ax10[1].set_ylim(-2e-4,2e-4) #not for 20 & 40um
        
        pps+=1
    else:
        continue


print(fnameall1)
#to heavy for pdf format
fig3.savefig(path+'orbit_distribution_correct_it1_{0}seeds.png'.format(eseed))
#fig5.savefig(path+'orbit_distribution_correct_it2_{0}seeds.png'.format(eseed))
#fig9.savefig(path+'orbit_distribution_tune_match_first_it_{0}seeds.png'.format(eseed))
#fig10.savefig(path+'orbit_distribution_tune_match_last_it_{0}seeds.png'.format(eseed))
    
max_orbit_x=np.max(orbit_x_all_seed)
max_orbit_y=np.max(orbit_y_all_seed)

print('analytical rms: max_orbit_x = ',max_orbit_x)
print('analytical rms: max_orbit_y = ',max_orbit_y)
print('\n')

print('rms_dist_x = ',rms_dist_x)
print('rms_dist_y = ',rms_dist_y)
print('\n')
"""
print('numerical rms x = ',rms_dist_x[0])
print('numerical rms y = ',rms_dist_y[0])
print('\n')
"""
mean_rms_dist_x2=np.mean(rms_dist_x)
mean_rms_dist_y2=np.mean(rms_dist_y)

print('numerical rms: mean_rms_dist_x = ',mean_rms_dist_x2)
print('numerical rms: mean_rms_dist_y = ',mean_rms_dist_y2)

#to heavy for pdf format
fig3.savefig(path+'orbit_distribution_correct_it1_{0}seeds.png'.format(eseed))
plt.show()
sys.exit()
#fig5.savefig(path+'orbit_distribution_correct_it2_{0}seeds.png'.format(eseed))
#fig9.savefig(path+'orbit_distribution_tune_match_first_it_{0}seeds.png'.format(eseed))
#fig10.savefig(path+'orbit_distribution_tune_match_last_it_{0}seeds.png'.format(eseed))

#plot of the rms: it without tune match
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x, ".", label='it n°1')
#ax2[0].plot(xseed2, rms_dist_x2, ".", label='it n°2')
ax2[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax2[0].set_ylabel("rms$_x$ [m]")
#ax2[0].set_ylim(0,500e-6) #for 80um
#ax2[0].set_ylim(0,6e-6) #for 20um
#ax2[0].set_ylim(0,1.9e-5) #for 40um
ax2[0].set_ylim(0,250e-6) #for 60um
#ax2[0].set_ylim(0,max_orbit_x+15e-5)
ax2[0].legend(fontsize=10,loc='best')
ax2[1].plot(xseed, rms_dist_y, ".", label='it n°1')
#ax2[1].plot(xseed2, rms_dist_y2, ".", label='it n°2')
ax2[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m]")
#ax2[1].set_ylim(0,500e-6) #for 80um
#ax2[1].set_ylim(4e-6,6e-6) #for 20um
#ax2[1].set_ylim(0.9e-5,1.9e-5) #for 40um
ax2[1].set_ylim(0,250e-6) #for 60um
#ax2[1].set_ylim(0,max_orbit_y+15e-5)
ax2[1].legend(fontsize=10,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_orbit_correct_{0}seeds.pdf'.format(eseed))

#plot of the mean: it without tune match
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x, ".", label='it n°1')
#ax4[0].plot(xseed2, ave_dist_x2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [m]")
#ax4[0].set_ylim(-50e-6,50e-6) #for 80um
#ax4[0].set_ylim(-3e-7,3e-7) #for 20um
#ax4[0].set_ylim(-6e-7,6e-7) #for 40um
ax4[0].set_ylim(-15e-7,15e-7) #for 60um
ax4[0].legend(fontsize=10,loc='best')
ax4[1].plot(xseed, ave_dist_y, ".", label='it n°1')
#ax4[1].plot(xseed2, ave_dist_y2, ".", label='it n°2')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m]")
#ax4[1].set_ylim(-50e-6,50e-6) #for 80um
#ax4[1].set_ylim(-3e-7,3e-7) #for 20um
#ax4[1].set_ylim(-6e-7,6e-7) #for 40um
ax4[1].set_ylim(-15e-7,15e-7) #for 60um
ax4[1].legend(fontsize=10,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_orbit_correct_{0}seeds.pdf'.format(eseed))
sys.exit()
#plot of the rms: it with tune match
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(xseed3, rms_dist_x3, ".", label='first it')
ax12[0].plot(xseed4, rms_dist_x4, ".", label='last it')
ax12[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax12[0].set_ylabel("rms$_x$ [m]")
#ax12[0].set_ylim(0,500e-6) #for 80um
#ax12[0].set_ylim(0,6e-6) #for 20um
#ax12[0].set_ylim(0,1.9e-5) #for 40um
ax12[0].set_ylim(0,250e-6) #for 60um
#ax12[0].set_ylim(0,max_orbit_x+15e-5)
ax12[0].legend(fontsize=10,loc='best')
ax12[1].plot(xseed3, rms_dist_y3, ".", label='first it')
ax12[1].plot(xseed4, rms_dist_y4, ".", label='last it')
ax12[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax12[1].set_xlabel("seed")
ax12[1].set_ylabel("rms$_y$ [m]")
#ax12[1].set_ylim(0,500e-6) #for 80um
#ax12[1].set_ylim(4e-6,6e-6) #for 20um
#ax12[1].set_ylim(0.9e-5,1.9e-5) #for 40um
ax12[1].set_ylim(0,250e-6) #for 60um
#ax12[1].set_ylim(0,max_orbit_y+15e-5)
ax12[1].legend(fontsize=10,loc='best')
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig12.savefig(path+'rms_orbit_tune_match_{0}seeds.pdf'.format(eseed))

#plot of the mean: it with tune match
fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax13[0].plot(xseed3, ave_dist_x3, ".", label='first it')
ax13[0].plot(xseed4, ave_dist_x4, ".", label='last it')
ax13[0].set_ylabel("mean$_x$ [m]")
#ax13[0].set_ylim(-50e-6,50e-6) #for 80um
#ax13[0].set_ylim(-3e-7,3e-7) #for 20um
#ax13[0].set_ylim(-6e-7,6e-7) #for 40um
ax13[0].set_ylim(-15e-7,15e-7) #for 60um
ax13[0].legend(fontsize=10,loc='best')
ax13[1].plot(xseed3, ave_dist_y3, ".", label='first it')
ax13[1].plot(xseed4, ave_dist_y4, ".", label='last it')
ax13[1].set_xlabel("seed")
ax13[1].set_ylabel("mean$_y$ [m]")
#ax13[1].set_ylim(-50e-6,50e-6) #for 80um
#ax13[1].set_ylim(-3e-7,3e-7) #for 20um
#ax13[1].set_ylim(-6e-7,6e-7) #for 40um
ax13[1].set_ylim(-15e-7,15e-7) #for 60um
ax13[1].legend(fontsize=10,loc='best')
fig13, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig13.savefig(path+'mean_orbit_tune_match_{0}seeds.pdf'.format(eseed))

#for iteration n°1: it without tune match
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
fig1.savefig(path+'hist_orbit_correct_it1_{0}seeds.pdf'.format(eseed))

#for iteration n°2: it without tune match
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
fig8.savefig(path+'hist_orbit_correct_it2_{0}seeds.pdf'.format(eseed))

#with tune match: first it
fig14, ax14=plt.subplots(nrows=2, ncols=2, sharey=True)
ax14[0,0].hist(rms_dist_x3,100)
ax14[0,0].set_xlabel("rms$_x$ [m]")
ax14[0,1].hist(ave_dist_x3,100)
ax14[0,1].set_xlabel("mean$_x$ [m]")
ax14[0,0].set_ylabel("counts")
ax14[1,0].hist(rms_dist_y3,100)
ax14[1,0].set_xlabel("rms$_y$ [m]")
ax14[1,1].hist(ave_dist_y3,100)
ax14[1,1].set_xlabel("mean$_y$ [m]")
ax14[1,0].set_ylabel("counts")
fig14, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig14.tight_layout()
fig14.savefig(path+'hist_orbit_tune_match_first_it_{0}seeds.pdf'.format(eseed))

#with tune match: last it
fig15, ax15=plt.subplots(nrows=2, ncols=2, sharey=True)
ax15[0,0].hist(rms_dist_x4,100)
ax15[0,0].set_xlabel("rms$_x$ [m]")
ax15[0,1].hist(ave_dist_x4,100)
ax15[0,1].set_xlabel("mean$_x$ [m]")
ax15[0,0].set_ylabel("counts")
ax15[1,0].hist(rms_dist_y4,100)
ax15[1,0].set_xlabel("rms$_y$ [m]")
ax15[1,1].hist(ave_dist_y4,100)
ax15[1,1].set_xlabel("mean$_y$ [m]")
ax15[1,0].set_ylabel("counts")
fig15, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig15.tight_layout()
fig15.savefig(path+'hist_orbit_tune_match_last_it_{0}seeds.pdf'.format(eseed))




