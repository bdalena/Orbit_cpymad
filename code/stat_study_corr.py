#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the correctors strength after correction and tune match.
'''
import math
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os

err_mq=150 #quads offset
err_ms=200 #sextupoles offset
err_bpm=10 #bpm offset
eseed=100 #nb of seeds

#path definition
#path="/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_100seeds_corrhplus/".format(err_mq)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_dxdy_{0}_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_ms)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_ms_dxdy_150_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_b1r_dpsi_mq_dxdy_{0}_tm_100seeds/'.format(err_mq)

fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig3, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig9, ax9=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig9, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig10, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig11, ax11=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig11, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig17, ax17=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig17, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig18, ax18=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig18, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

#line sextuoff it1
rms_dist_x_it1=np.empty(0)
ave_dist_x_it1=np.empty(0)
rms_dist_y_it1=np.empty(0)
ave_dist_y_it1=np.empty(0)

#line sextuoff it2
rms_dist_x_it2=np.empty(0)
ave_dist_x_it2=np.empty(0)
rms_dist_y_it2=np.empty(0)
ave_dist_y_it2=np.empty(0)

#ring sextuoff last it
rms_dist_x_it4=np.empty(0)
ave_dist_x_it4=np.empty(0)
rms_dist_y_it4=np.empty(0)
ave_dist_y_it4=np.empty(0)

#ring sextuon
rms_dist_x_it5=np.empty(0)
ave_dist_x_it5=np.empty(0)
rms_dist_y_it5=np.empty(0)
ave_dist_y_it5=np.empty(0)

#tm first it
rms_dist_x_it6=np.empty(0)
ave_dist_x_it6=np.empty(0)
rms_dist_y_it6=np.empty(0)
ave_dist_y_it6=np.empty(0)

#tm last it
rms_dist_x_it7=np.empty(0)
ave_dist_x_it7=np.empty(0)
rms_dist_y_it7=np.empty(0)
ave_dist_y_it7=np.empty(0)

xseed=np.empty(0)
xseed3 = np.empty(0)
xseed4 = np.empty(0)

iis=0
kks=0
pps=0

p=182 #Gev
Brho=3.3356*p

mean_corr_x_all_seed=np.empty(0)
mean_corr_y_all_seed=np.empty(0)

for j in range(eseed):
    
    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(j+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(j+1)
    mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)
    mean_corr_x_all_seed=np.append(mean_corr_x_all_seed,mean_corr_x)
    mean_corr_y_all_seed=np.append(mean_corr_y_all_seed,mean_corr_y)

max_corr_x3=3*np.max(mean_corr_x_all_seed)*Brho
max_corr_y3=3*np.max(mean_corr_y_all_seed)*Brho
max_corr_x=np.max(mean_corr_x_all_seed)*Brho
max_corr_y=np.max(mean_corr_y_all_seed)*Brho

print('analytical rms: max_corr_x = ',max_corr_x)
print('analytical rms: max_corr_y = ',max_corr_y)
print('\n')

for i in range(eseed):

    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(i+1)

    #if os.path.exists(path+"test_seed{0}/cx_fccee_heb_mic_all_line_sextuoff_it2.tab".format(i+1)):
    if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs".format(i+1)):
        fnameallx=path+"test_seed{0}/cx_fccee_heb_mic_all_line_sextuoff_it2.tab".format(i+1)
        fnameally=path+"test_seed{0}/cy_fccee_heb_mic_all_line_sextuoff_it2.tab".format(i+1)
    
    if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs".format(i+1)):
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
        if i==1:
            ax3[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax3[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax3[0].legend(fontsize=8,loc='best')
        ax3[0].set_ylim(-300e-4,300e-4) #for 200 um
        #ax3[0].set_ylim(-180e-4,180e-4) #for 100 um
        #ax3[0].set_ylim(-120e-4,120e-4) #for 80 um
        ax3[1].plot(ycorr, optics_ally['PY.OLD']*Brho, ".")
        ax3[1].set_xlabel("Correctors")
        ax3[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax3[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax3[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax3[1].legend(fontsize=8,loc='best')
        #ax3[1].set_ylim(-120e-4,120e-4)  #for 80 um
        #ax3[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax3[1].set_ylim(-300e-4,300e-4) #for 200 um
        
        ax5[0].plot(xcorr, optics_allx['PX.CORRECTION']*Brho, ".")
        ax5[0].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax5[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax5[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax5[0].legend(fontsize=8,loc='best')
        #ax5[0].set_ylim(-180e-4,180e-4) #for 100 um
        ax5[0].set_ylim(-300e-4,300e-4) #for 200 um
        ax5[1].plot(ycorr, optics_ally['PY.CORRECTION']*Brho, ".")
        ax5[1].set_xlabel("Correctors")
        ax5[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax5[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax5[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax5[1].legend(fontsize=8,loc='best')
        #ax5[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax5[1].set_ylim(-300e-4,300e-4) #for 200 um
        
        iis+=1
    else:
        continue

    if os.path.exists(path+"test_seed{0}/cx_fccee_heb_mic_all_ring_sextuon_it1.tab".format(i+1)):
        fnameallx3=path+"test_seed{0}/cx_fccee_heb_mic_all_ring_sextuon_it1.tab".format(i+1)
        fnameally3=path+"test_seed{0}/cy_fccee_heb_mic_all_ring_sextuon_it1.tab".format(i+1)

    if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(i+1)):
        xseed3=np.append(xseed3,kks)
        head_opt=pd.read_csv(fnameallx3, header=6, sep='\s+', nrows=0).columns[1:]
        optics_allx3=pd.read_csv(fnameallx3, skiprows=8, sep='\s+', names=head_opt)
        optics_allx3=optics_allx3.reset_index(drop=True)
        xcorr3=np.arange(0,len(optics_allx3),1)
        rms_dist_x_it4=np.append(rms_dist_x_it4,np.std(optics_allx3['PX.OLD']*Brho))
        ave_dist_x_it4=np.append(ave_dist_x_it4,np.mean(optics_allx3['PX.OLD']*Brho))
        rms_dist_x_it5=np.append(rms_dist_x_it5,np.std(optics_allx3['PX.CORRECTION']*Brho))
        ave_dist_x_it5=np.append(ave_dist_x_it5,np.mean(optics_allx3['PX.CORRECTION']*Brho))

        head_opt=pd.read_csv(fnameally3, header=6, sep='\s+', nrows=0).columns[1:]
        optics_ally3=pd.read_csv(fnameally3, skiprows=8, sep='\s+', names=head_opt)
        optics_ally3=optics_ally3.reset_index(drop=True)
        ycorr3=np.arange(0,len(optics_ally3),1)
        rms_dist_y_it4=np.append(rms_dist_y_it4,np.std(optics_ally3['PY.OLD']*Brho))
        ave_dist_y_it4=np.append(ave_dist_y_it4,np.mean(optics_ally3['PY.OLD']*Brho))
        rms_dist_y_it5=np.append(rms_dist_y_it5,np.std(optics_ally3['PY.CORRECTION']*Brho))
        ave_dist_y_it5=np.append(ave_dist_y_it5,np.mean(optics_ally3['PY.CORRECTION']*Brho))
        
        ax10[0].plot(xcorr3, optics_allx3['PX.OLD']*Brho, ".")
        ax10[0].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax10[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax10[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax10[0].legend(fontsize=8,loc='best')
        #ax10[0].set_ylim(-180e-4,180e-4) #for 100 um
        ax10[0].set_ylim(-250e-4,250e-4) #for 200 um 
        ax10[1].plot(ycorr3, optics_ally3['PY.OLD']*Brho, ".")
        ax10[1].set_xlabel("Correctors")
        ax10[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax10[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax10[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax10[1].legend(fontsize=8,loc='best')
         #ax10[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax10[1].set_ylim(-250e-4,250e-4) #for 200 um

        ax11[0].plot(xcorr3, optics_allx3['PX.CORRECTION']*Brho, ".")
        ax11[0].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax11[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax11[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax11[0].legend(fontsize=8,loc='best')
        #ax11[0].set_ylim(-180e-4,180e-4) #for 100 um
        ax11[0].set_ylim(-250e-4,250e-4) #for 200 um
        ax11[1].plot(ycorr3, optics_ally3['PY.CORRECTION']*Brho, ".")
        ax11[1].set_xlabel("Correctors")
        ax11[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax11[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax11[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax11[1].legend(fontsize=8,loc='best')
        #ax11[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax11[1].set_ylim(-250e-4,250e-4) #for 200 um
        
        kks+=1
    else:
        continue

    if os.path.exists(path+"test_seed{0}/cx_fccee_heb_mic_all_tm_it1.tab".format(i+1)):
        fnameallx4=path+"test_seed{0}/cx_fccee_heb_mic_all_tm_it1.tab".format(i+1)
        fnameally4=path+"test_seed{0}/cy_fccee_heb_mic_all_tm_it1.tab".format(i+1)

    if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(i+1)):
        xseed4=np.append(xseed4,pps)
        head_opt=pd.read_csv(fnameallx4, header=6, sep='\s+', nrows=0).columns[1:]
        optics_allx4=pd.read_csv(fnameallx4, skiprows=8, sep='\s+', names=head_opt)
        optics_allx4=optics_allx4.reset_index(drop=True)
        xcorr4=np.arange(0,len(optics_allx4),1)
        rms_dist_x_it6=np.append(rms_dist_x_it6,np.std(optics_allx4['PX.OLD']*Brho))
        ave_dist_x_it6=np.append(ave_dist_x_it6,np.mean(optics_allx4['PX.OLD']*Brho))
        rms_dist_x_it7=np.append(rms_dist_x_it7,np.std(optics_allx4['PX.CORRECTION']*Brho))
        ave_dist_x_it7=np.append(ave_dist_x_it7,np.mean(optics_allx4['PX.CORRECTION']*Brho))

        head_opt=pd.read_csv(fnameally4, header=6, sep='\s+', nrows=0).columns[1:]
        optics_ally4=pd.read_csv(fnameally4, skiprows=8, sep='\s+', names=head_opt)
        optics_ally4=optics_ally4.reset_index(drop=True)
        ycorr4=np.arange(0,len(optics_ally4),1)
        rms_dist_y_it6=np.append(rms_dist_y_it6,np.std(optics_ally4['PY.OLD']*Brho))
        ave_dist_y_it6=np.append(ave_dist_y_it6,np.mean(optics_ally4['PY.OLD']*Brho))
        rms_dist_y_it7=np.append(rms_dist_y_it7,np.std(optics_ally4['PY.CORRECTION']*Brho))
        ave_dist_y_it7=np.append(ave_dist_y_it7,np.mean(optics_ally4['PY.CORRECTION']*Brho))
        
        ax17[0].plot(xcorr4, optics_allx4['PX.OLD']*Brho, ".")
        ax17[0].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax17[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax17[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax17[0].legend(fontsize=8,loc='best')
        #ax17[0].set_ylim(-180e-4,180e-4) #for 100 um
        ax17[0].set_ylim(-250e-4,250e-4) #for 200 um 
        ax17[1].plot(ycorr4, optics_ally4['PY.OLD']*Brho, ".")
        ax17[1].set_xlabel("Correctors")
        ax17[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax17[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax17[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax17[1].legend(fontsize=8,loc='best')
        #ax17[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax17[1].set_ylim(-250e-4,250e-4) #for 200 um

        ax18[0].plot(xcorr4, optics_allx4['PX.CORRECTION']*Brho, ".")
        ax18[0].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax18[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='3 analytical rms')
            ax18[0].axhline(y=-max_corr_x3, color='r', linestyle='--')
            ax18[0].legend(fontsize=8,loc='best')
        #ax18[0].set_ylim(-180e-4,180e-4) #for 100 um
        ax18[0].set_ylim(-250e-4,250e-4) #for 200 um
        ax18[1].plot(ycorr4, optics_ally4['PY.CORRECTION']*Brho, ".")
        ax18[1].set_xlabel("Correctors")
        ax18[1].set_ylabel("Corrector strength [Tm]")
        if i==1:
            ax18[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='3 analytical rms')
            ax18[1].axhline(y=-max_corr_y3, color='r', linestyle='--')
            ax18[1].legend(fontsize=8,loc='best')
        #ax18[1].set_ylim(-180e-4,180e-4) #for 100 um
        ax18[1].set_ylim(-250e-4,250e-4) #for 200 um
        
        pps+=1
    else:
        continue
    
print('numerical rms x = ',rms_dist_x_it2[0])
print('numerical rms y = ',rms_dist_y_it2[0])
print('\n')

mean_rms_dist_x2=np.mean(rms_dist_x_it2)
mean_rms_dist_y2=np.mean(rms_dist_y_it2)

print('numerical rms: mean_rms_dist_x2 = ',mean_rms_dist_x2)
print('numerical rms: mean_rms_dist_y2 = ',mean_rms_dist_y2)

sys.exit()

#to heavy for pdf format 
fig3.savefig(path+'correctors_strength_distribution_line_sextuoff_it1_{0}seeds.png'.format(eseed))
fig5.savefig(path+'correctors_strength_distribution_line_sextuoff_it2_{0}seeds.png'.format(eseed))   
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):    
    fig10.savefig(path+'correctors_strength_distribution_ring_sextuoff_last_it_{0}seeds.png'.format(eseed))
    fig11.savefig(path+'correctors_strength_distribution_ring_sextuon_it1_{0}seeds.png'.format(eseed))
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(eseed)):    
    fig17.savefig(path+'correctors_strength_distribution_tm_it0_{0}seeds.png'.format(eseed))
    fig18.savefig(path+'correctors_strength_distribution_tm_it1_{0}seeds.png'.format(eseed))

#plot of the rms: sextuoff line
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x_it1, ".", label='it n°1')
ax2[0].plot(xseed, rms_dist_x_it2, ".", label='it n°2')
#ax2[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='analytical 3 rms')
ax2[0].axhline(y=max_corr_x, color='g', linestyle='--', label='analytical rms')
ax2[0].set_ylabel("rms$_x$ [Tm]")
ax2[0].set_ylim(0.003,max_corr_x+80e-5)
ax2[0].legend(fontsize=8,loc='best')
ax2[1].plot(xseed, rms_dist_y_it1, ".", label='it n°1')
ax2[1].plot(xseed, rms_dist_y_it2, ".", label='it n°2')
#ax2[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='analytical 3 rms')
ax2[1].axhline(y=max_corr_y, color='g', linestyle='--', label='analytical rms')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [Tm]")
ax2[1].set_ylim(0.003,max_corr_y+80e-5)
ax2[1].legend(fontsize=8,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_correctors_strength_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the mean: sextuoff line
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x_it1, ".", label='it n°1')
ax4[0].plot(xseed, ave_dist_x_it2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [Tm]")
ax4[0].set_ylim(-200e-6,200e-6)
ax4[0].legend(fontsize=8,loc='best')
ax4[1].plot(xseed, ave_dist_y_it1, ".", label='it n°1')
ax4[1].plot(xseed, ave_dist_y_it2, ".", label='it n°1')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [Tm]")
ax4[1].set_ylim(-200e-6,200e-6)
ax4[1].legend(fontsize=8,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_correctors_strength_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the rms: sextuoff/on ring
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):
    fig11, ax11=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    #ax11[0].plot(xseed2, rms_dist_x_it3, ".", label='MS_off: first it')
    ax11[0].plot(xseed3, rms_dist_x_it4, ".", label='MS_off: last it')
    ax11[0].plot(xseed3, rms_dist_x_it5, ".", label='MS_on')
    #ax11[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='analytical 3 rms')
    ax11[0].axhline(y=max_corr_x, color='g', linestyle='--', label='analytical rms')
    ax11[0].set_ylabel("rms$_x$ [Tm]")
    ax11[0].set_ylim(0.003,max_corr_x+80e-5)
    ax11[0].legend(fontsize=8,loc='best')
    #ax11[1].plot(xseed2, rms_dist_y_it3, ".", label='MS_off: first it')
    ax11[1].plot(xseed3, rms_dist_y_it4, ".", label='MS_off: last it')
    ax11[1].plot(xseed3, rms_dist_y_it5, ".", label='MS_on')
    #ax11[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='analytical 3 rms')
    ax11[1].axhline(y=max_corr_y, color='g', linestyle='--', label='analytical 1 rms')
    ax11[1].set_xlabel("seed")
    ax11[1].set_ylabel("rms$_y$ [Tm]")
    ax11[1].set_ylim(0.003,max_corr_y+80e-5)
    ax11[1].legend(fontsize=8,loc='best')
    fig11, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
    fig11.savefig(path+'rms_correctors_strength_ring_{0}seeds.pdf'.format(eseed))    

#plot of the mean: sextuoff/on ring
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):
    fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax12[0].plot(xseed3, ave_dist_x_it4, ".", label='MS_off: last it')
    ax12[0].plot(xseed3, ave_dist_x_it5, ".", label='MS_on')
    ax12[0].set_ylabel("mean$_x$ [Tm]")
    ax12[0].set_ylim(-200e-6,200e-6)
    ax12[0].legend(fontsize=8,loc='best')
    ax12[1].plot(xseed3, ave_dist_y_it4, ".", label='MS_off: last it')
    ax12[1].plot(xseed3, ave_dist_y_it5, ".", label='MS_on')
    ax12[1].set_xlabel("seed")
    ax12[1].set_ylabel("mean$_y$ [Tm]")
    ax12[1].set_ylim(-200e-6,200e-6)
    ax12[1].legend(fontsize=8,loc='best')
    fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
    fig12.savefig(path+'mean_correctors_strength_ring_{0}seeds.pdf'.format(eseed))

#plot of the rms: tm 
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(eseed)):
    fig20, ax20=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax20[0].plot(xseed4, rms_dist_x_it6, ".", label='tm first it')
    ax20[0].plot(xseed4, rms_dist_x_it7, ".", label='tm last it')
    #ax20[0].axhline(y=max_corr_x3, color='r', linestyle='--', label='analytical 3 rms')
    ax20[0].axhline(y=max_corr_x, color='g', linestyle='--', label='analytical rms')
    ax20[0].set_ylabel("rms$_x$ [Tm]")
    ax20[0].set_ylim(0.003,max_corr_x+80e-5)
    ax20[0].legend(fontsize=8,loc='best')
    ax20[1].plot(xseed4, rms_dist_y_it6, ".", label='tm first it')
    ax20[1].plot(xseed4, rms_dist_y_it7, ".", label='tm last it')
    #ax20[1].axhline(y=max_corr_y3, color='r', linestyle='--', label='analytical 3 rms')
    ax20[1].axhline(y=max_corr_y, color='g', linestyle='--', label='analytical 1 rms')
    ax20[1].set_xlabel("seed")
    ax20[1].set_ylabel("rms$_y$ [Tm]")
    ax20[1].set_ylim(0.003,max_corr_y+80e-5)
    ax20[1].legend(fontsize=8,loc='best')
    fig20, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
    fig20.savefig(path+'rms_correctors_strength_tm_{0}seeds.pdf'.format(eseed))   

#plot of the mean: tm 
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(eseed)):
    fig19, ax19=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax19[0].plot(xseed4, ave_dist_x_it6, ".", label='tm first it')
    ax19[0].plot(xseed4, ave_dist_x_it7, ".", label='tm last it')
    ax19[0].set_ylabel("mean$_x$ [Tm]")
    ax19[0].set_ylim(-200e-6,200e-6)
    ax19[0].legend(fontsize=8,loc='best')
    ax19[1].plot(xseed4, ave_dist_y_it6, ".", label='tm first it')
    ax19[1].plot(xseed4, ave_dist_y_it7, ".", label='tm last it')
    ax19[1].set_xlabel("seed")
    ax19[1].set_ylabel("mean$_y$ [Tm]")
    ax19[1].set_ylim(-200e-6,200e-6)
    ax19[1].legend(fontsize=8,loc='best')
    fig19, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
    fig19.savefig(path+'mean_correctors_strength_tm_{0}seeds.pdf'.format(eseed))

    
'''
#line sextuoff it1
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
fig1.savefig(path+'hist_correctors_strength_line_sextuoff_it1_{0}seeds.pdf'.format(eseed))

#line sextuoff it2
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
fig8.savefig(path+'hist_correctors_strength_line_sextuoff_it2_{0}seeds.pdf'.format(eseed))

#ring sextuoff first it
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs".format(eseed)):
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
    fig13.savefig(path+'hist_correctors_strength_ring_sextuoff_first_it_{0}seeds.pdf'.format(eseed))

#ring sextuoff last it
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):
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
    fig14.savefig(path+'hist_correctors_strength_ring_sextuoff_last_it_{0}seeds.pdf'.format(eseed))

#ring sextuon
if os.path.exists(path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):
    fig16, ax16=plt.subplots(nrows=2, ncols=2, sharey=True)
    ax16[0,0].hist(rms_dist_x_it5,100)
    ax16[0,0].set_xlabel("rms$_x$ [Tm]")
    ax16[0,1].hist(ave_dist_x_it5,100)
    ax16[0,1].set_xlabel("mean$_x$ [Tm]")
    ax16[0,0].set_ylabel("counts")
    ax16[1,0].hist(rms_dist_y_it5,100)
    ax16[1,0].set_xlabel("rms$_y$ [Tm]")
    ax16[1,1].hist(ave_dist_y_it5,100)
    ax16[1,1].set_xlabel("mean$_y$ [Tm]")
    ax16[1,0].set_ylabel("counts")
    fig16, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
    fig16.tight_layout()
    fig16.savefig(path+'hist_orbit_ring_sextuon_it1_{0}seeds.pdf'.format(eseed))
'''
