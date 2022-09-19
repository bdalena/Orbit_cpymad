#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the correctors strength after correction and tune match.
'''

import matplotlib.pyplot as plt
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os

err_mq=90 #quads offset
eseed=100 #nb of seeds

#path definition
#path='/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_100seeeds_corrhplus/'.format(err_mq)
#path='/home/td271008/work/cpymadx/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_2it_100seeeds_corrhplus/'.format(err_mq)
path='/home/td271008/work/cpymadx/Orbit_cpymad/tune_match_mb_fielderr_roll_mq_offset_{0}_IP5_100seeeds_corrhplus/'.format(err_mq)
#path='/home/td271008/work/cpymadx/Orbit_cpymad/tune_match_mb_fielderr_roll_mq_offset_{0}_IP5_10seeeds_corrhplus/'.format(err_mq)

fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig9, ax9=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig9, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig10, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

#CORRECT: it 1
rms_dist_x_it1=np.empty(0)
ave_dist_x_it1=np.empty(0)
rms_dist_y_it1=np.empty(0)
ave_dist_y_it1=np.empty(0)

#CORRECT: it 2
rms_dist_x_it2=np.empty(0)
ave_dist_x_it2=np.empty(0)
rms_dist_y_it2=np.empty(0)
ave_dist_y_it2=np.empty(0)

#tune match: first it
rms_dist_x_it3=np.empty(0)
ave_dist_x_it3=np.empty(0)
rms_dist_y_it3=np.empty(0)
ave_dist_y_it3=np.empty(0)

#tune match: last it
rms_dist_x_it4=np.empty(0)
ave_dist_x_it4=np.empty(0)
rms_dist_y_it4=np.empty(0)
ave_dist_y_it4=np.empty(0)

xseed=np.empty(0)
xseed2=np.empty(0)
xseed3 = np.empty(0)
xseed4 = np.empty(0)

iis=0
jjs=0
kks=0
pps=0

p=182 #Gev
Brho=3.3356*p

mean_corr_x_all_seed=np.empty(0)
mean_corr_y_all_seed=np.empty(0)

for i in range(eseed):

    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(i+1)

    mean_corr_x_seed,mean_corr_y_seed,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)
    mean_corr_x_all_seed=np.append(mean_corr_x_all_seed,mean_corr_x_seed)
    mean_corr_y_all_seed=np.append(mean_corr_y_all_seed,mean_corr_y_seed)
    
    fnameallx=path+"test_seed{0}/cx_fccee_heb_mic_all_line_it2.tab".format(i+1)
    fnameally=path+"test_seed{0}/cy_fccee_heb_mic_all_line_it2.tab".format(i+1)
    
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
        ax3[0].set_ylim(-180e-4,180e-4) #for 100 um
        #ax3[0].set_ylim(-120e-4,120e-4) #for 80 um
        ax3[1].plot(ycorr, optics_ally['PY.OLD']*Brho, ".")
        ax3[1].set_xlabel("Correctors")
        ax3[1].set_ylabel("Corrector strength [Tm]")
        #ax3[1].set_ylim(-120e-4,120e-4)  #for 80 um
        ax3[1].set_ylim(-180e-4,180e-4) #for 100 um
        
        ax5[0].plot(xcorr, optics_allx['PX.CORRECTION']*Brho, ".")
        ax5[0].set_ylabel("Corrector strength [Tm]")
        #ax5[0].set_ylim(-120e-4,120e-4)  #for 80 um
        ax5[0].set_ylim(-180e-4,180e-4) #for 100 um 
        ax5[1].plot(ycorr, optics_ally['PY.CORRECTION']*Brho, ".")
        ax5[1].set_xlabel("Correctors")
        ax5[1].set_ylabel("Corrector strength [Tm]")
        #ax5[1].set_ylim(-120e-4,120e-4)  #for 80 um
        ax5[1].set_ylim(-180e-4,180e-4) #for 100 um 
        
        iis+=1

    fnameallx2=path+"test_seed{0}/cx_fccee_heb_mic_all_ring_it1.tab".format(i+1)
    fnameally2=path+"test_seed{0}/cy_fccee_heb_mic_all_ring_it1.tab".format(i+1)

    if os.path.exists(fnameallx2):
        xseed2=np.append(xseed2,jjs)
        head_opt=pd.read_csv(fnameallx2, header=6, sep='\s+', nrows=0).columns[1:]
        optics_allx2=pd.read_csv(fnameallx2, skiprows=8, sep='\s+', names=head_opt)
        optics_allx2=optics_allx2.reset_index(drop=True)
        xcorr2=np.arange(0,len(optics_allx2),1)
        rms_dist_x_it3=np.append(rms_dist_x_it3,np.std(optics_allx2['PX.CORRECTION']*Brho))
        ave_dist_x_it3=np.append(ave_dist_x_it3,np.mean(optics_allx2['PX.CORRECTION']*Brho))

        head_opt=pd.read_csv(fnameally2, header=6, sep='\s+', nrows=0).columns[1:]
        optics_ally2=pd.read_csv(fnameally2, skiprows=8, sep='\s+', names=head_opt)
        optics_ally2=optics_ally2.reset_index(drop=True)
        ycorr2=np.arange(0,len(optics_ally2),1)
        rms_dist_y_it3=np.append(rms_dist_y_it3,np.std(optics_ally2['PY.CORRECTION']*Brho))
        ave_dist_y_it3=np.append(ave_dist_y_it3,np.mean(optics_ally2['PY.CORRECTION']*Brho))
        
        ax9[0].plot(xcorr2, optics_allx2['PX.CORRECTION']*Brho, ".")
        ax9[0].set_ylabel("Corrector strength [Tm]")
        #ax9[0].set_ylim(-120e-4,120e-4)  #for 80 um
        ax9[0].set_ylim(-180e-4,180e-4) #for 100 um 
        ax9[1].plot(ycorr2, optics_ally2['PY.CORRECTION']*Brho, ".")
        ax9[1].set_xlabel("Correctors")
        ax9[1].set_ylabel("Corrector strength [Tm]")
        #ax9[1].set_ylim(-120e-4,120e-4)  #for 80 um
        ax9[1].set_ylim(-180e-4,180e-4) #for 100 um 
        
        jjs+=1

    cnt=1
    fnameallx3=path+"test_seed{0}/cx_fccee_heb_mic_all_ring_it{1}.tab".format(i+1,cnt)
    
    while os.path.exists(fnameallx3):
        fnameallx3=path+"test_seed{0}/cx_fccee_heb_mic_all_ring_it{1}.tab".format(i+1,cnt+1)
        cnt+=1
    
    fnameallx3=path+"test_seed{0}/cx_fccee_heb_mic_all_ring_it{1}.tab".format(i+1,cnt-1)
    fnameally3=path+"test_seed{0}/cy_fccee_heb_mic_all_ring_it{1}.tab".format(i+1,cnt-1)

    if os.path.exists(fnameallx3):
        xseed3=np.append(xseed3,kks)
        head_opt=pd.read_csv(fnameallx3, header=6, sep='\s+', nrows=0).columns[1:]
        optics_allx3=pd.read_csv(fnameallx3, skiprows=8, sep='\s+', names=head_opt)
        optics_allx3=optics_allx3.reset_index(drop=True)
        xcorr3=np.arange(0,len(optics_allx3),1)
        rms_dist_x_it4=np.append(rms_dist_x_it4,np.std(optics_allx3['PX.CORRECTION']*Brho))
        ave_dist_x_it4=np.append(ave_dist_x_it4,np.mean(optics_allx3['PX.CORRECTION']*Brho))

        head_opt=pd.read_csv(fnameally3, header=6, sep='\s+', nrows=0).columns[1:]
        optics_ally3=pd.read_csv(fnameally2, skiprows=8, sep='\s+', names=head_opt)
        optics_ally3=optics_ally3.reset_index(drop=True)
        ycorr3=np.arange(0,len(optics_ally3),1)
        rms_dist_y_it4=np.append(rms_dist_y_it4,np.std(optics_ally3['PY.CORRECTION']*Brho))
        ave_dist_y_it4=np.append(ave_dist_y_it4,np.mean(optics_ally3['PY.CORRECTION']*Brho))
        
        ax10[0].plot(xcorr3, optics_allx3['PX.CORRECTION']*Brho, ".")
        ax10[0].set_ylabel("Corrector strength [Tm]")
        #ax10[0].set_ylim(-120e-4,120e-4)  #for 80 um
        ax10[0].set_ylim(-180e-4,180e-4) #for 100 um 
        ax10[1].plot(ycorr3, optics_ally3['PY.CORRECTION']*Brho, ".")
        ax10[1].set_xlabel("Correctors")
        ax10[1].set_ylabel("Corrector strength [Tm]")
        #ax10[1].set_ylim(-120e-4,120e-4)  #for 80 um
        ax10[1].set_ylim(-180e-4,180e-4) #for 100 um 
        
        kks+=1
    
max_corr_x3=3*np.max(mean_corr_x_all_seed)*Brho
max_corr_y3=3*np.max(mean_corr_y_all_seed)*Brho
max_corr_x=np.max(mean_corr_x_all_seed)*Brho
max_corr_y=np.max(mean_corr_y_all_seed)*Brho

#to heavy for pdf format 
fig3.savefig(path+'correctors_strength_distribution_it1_{0}seeds.png'.format(eseed))
fig5.savefig(path+'correctors_strength_distribution_it2_{0}seeds.png'.format(eseed))
fig9.savefig(path+'correctors_strength_distribution_tune_match_first_it_{0}seeds.png'.format(eseed))
fig10.savefig(path+'correctors_strength_distribution_tune_match_last_it_{0}seeds.png'.format(eseed))

#plot of the rms: it without tune match
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x_it1, ".", label='it n°1')
ax2[0].plot(xseed, rms_dist_x_it2, ".", label='it n°2')
ax2[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='analytical 3 rms')
ax2[0].axhline(y=max_corr_x, color='g', linestyle='--', label='analytical 1 rms')
ax2[0].set_ylabel("rms$_x$ [Tm]")
ax2[0].set_ylim(0,max_corr_x3+80e-5)
ax2[0].legend(fontsize=10,loc='best')
ax2[1].plot(xseed, rms_dist_y_it1, ".", label='it n°1')
ax2[1].plot(xseed, rms_dist_y_it2, ".", label='it n°2')
ax2[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='analytical 3 rms')
ax2[1].axhline(y=max_corr_y, color='g', linestyle='--', label='analytical 1 rms')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [Tm]")
ax2[1].set_ylim(0,max_corr_y3+80e-5)
ax2[1].legend(fontsize=10,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_correctors_strength_correct_{0}seeds.pdf'.format(eseed))

#plot of the mean: it without tune match
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x_it1, ".", label='it n°1')
ax4[0].plot(xseed, ave_dist_x_it2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [Tm]")
ax4[0].set_ylim(-300e-6,300e-6) #for 80, 90 & 100um
#ax4[0].set_ylim(-100e-6,100e-6) #for 20, 40 & 60um
ax4[0].legend(fontsize=10,loc='best')
ax4[1].plot(xseed, ave_dist_y_it1, ".", label='it n°1')
ax4[1].plot(xseed, ave_dist_y_it2, ".", label='it n°1')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [Tm]")
ax4[1].set_ylim(-300e-6,300e-6) #for 80, 90 & 100um
#ax4[1].set_ylim(-100e-6,100e-6) #for 20, 40 & 60um
ax4[1].legend(fontsize=10,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_correctors_strength_correct_{0}seeds.pdf'.format(eseed))

#plot of the rms: it with tune match
fig11, ax11=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax11[0].plot(xseed2, rms_dist_x_it3, ".", label='first it')
ax11[0].plot(xseed3, rms_dist_x_it4, ".", label='last it')
ax11[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='analytical 3 rms')
ax11[0].axhline(y=max_corr_x, color='g', linestyle='--', label='analytical 1 rms')
ax11[0].set_ylabel("rms$_x$ [Tm]")
ax11[0].set_ylim(0,max_corr_x3+80e-5)
ax11[0].legend(fontsize=10,loc='best')
ax11[1].plot(xseed2, rms_dist_y_it3, ".", label='first it')
ax11[1].plot(xseed3, rms_dist_y_it4, ".", label='last it')
ax11[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='analytical 3 rms')
ax11[1].axhline(y=max_corr_y, color='g', linestyle='--', label='analytical 1 rms')
ax11[1].set_xlabel("seed")
ax11[1].set_ylabel("rms$_y$ [Tm]")
ax11[1].set_ylim(0,max_corr_y3+80e-5)
ax11[1].legend(fontsize=10,loc='best')
fig11, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig11.savefig(path+'rms_correctors_strength_tune_match_{0}seeds.pdf'.format(eseed))

#plot of the mean: it with tune match
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(xseed2, ave_dist_x_it3, ".", label='first it')
ax12[0].plot(xseed3, ave_dist_x_it4, ".", label='last it')
ax12[0].set_ylabel("mean$_x$ [Tm]")
ax12[0].set_ylim(-300e-6,300e-6) #for 80, 90 & 100um
#ax12[0].set_ylim(-100e-6,100e-6) #for 20, 40 & 60um
ax12[0].legend(fontsize=10,loc='best')
ax12[1].plot(xseed2, ave_dist_y_it3, ".", label='first it')
ax12[1].plot(xseed3, ave_dist_y_it4, ".", label='last it')
ax12[1].set_xlabel("seed")
ax12[1].set_ylabel("mean$_y$ [Tm]")
ax12[1].set_ylim(-300e-6,300e-6) #for 80, 90 & 100um
#ax12[1].set_ylim(-100e-6,100e-6) #for 20, 40 & 60um
ax12[1].legend(fontsize=10,loc='best')
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig12.savefig(path+'mean_correctors_strength_tune_match_{0}seeds.pdf'.format(eseed))

#for iteration n°1: it without tune match
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
fig1.savefig(path+'hist_correctors_strength_it1_{0}seeds.pdf'.format(eseed))

#for iteration n°2: it without tune match
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
fig8.savefig(path+'hist_correctors_strength_it2_{0}seeds.pdf'.format(eseed))

#for iteration n°2: it with tune match
fig13, ax13=plt.subplots(nrows=2, ncols=2, sharey=True)
ax13[0,0].hist(rms_dist_x_it3,100)
ax13[0,0].set_xlabel("rms$_x$ [Tm]")
ax13[0,1].hist(ave_dist_x_it3,100)
ax13[0,1].set_xlabel("mean$_x$ [Tm]")
ax13[0,0].set_ylabel("counts")
ax13[1,0].hist(rms_dist_y_it3,100)
ax13[1,0].set_xlabel("rms$_y$ [Tm]")
ax13[1,1].hist(ave_dist_y_it3,100)
ax13[1,1].set_xlabel("mean$_y$ [Tm]")
ax13[1,0].set_ylabel("counts")
fig13, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig13.tight_layout()
fig13.savefig(path+'hist_correctors_strength_tune_match_first_it_{0}seeds.pdf'.format(eseed))

#for iteration n°3: it with tune match
fig14, ax14=plt.subplots(nrows=2, ncols=2, sharey=True)
ax14[0,0].hist(rms_dist_x_it4,100)
ax14[0,0].set_xlabel("rms$_x$ [Tm]")
ax14[0,1].hist(ave_dist_x_it4,100)
ax14[0,1].set_xlabel("mean$_x$ [Tm]")
ax14[0,0].set_ylabel("counts")
ax14[1,0].hist(rms_dist_y_it4,100)
ax14[1,0].set_xlabel("rms$_y$ [Tm]")
ax14[1,1].hist(ave_dist_y_it4,100)
ax14[1,1].set_xlabel("mean$_y$ [Tm]")
ax14[1,0].set_ylabel("counts")
fig14, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig14.tight_layout()
fig14.savefig(path+'hist_correctors_strength_tune_match_last_it_{0}seeds.pdf'.format(eseed))
