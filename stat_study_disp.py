#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the normalized dispersion after correction  and tune match.
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
rms_disp_x =np.empty(0)
ave_disp_x =np.empty(0)
rms_disp_y =np.empty(0)
ave_disp_y =np.empty(0)

#CORRECT: it 2
rms_disp_x2 =np.empty(0)
ave_disp_x2 =np.empty(0)
rms_disp_y2 =np.empty(0)
ave_disp_y2 =np.empty(0)

#tune match: first it
rms_disp_x3 =np.empty(0)
ave_disp_x3 =np.empty(0)
rms_disp_y3 =np.empty(0)
ave_disp_y3 =np.empty(0)

#tune match: last it
rms_disp_x4 =np.empty(0)
ave_disp_x4 =np.empty(0)
rms_disp_y4 =np.empty(0)
ave_disp_y4 =np.empty(0)

xseed = np.empty(0)
xseed2 = np.empty(0)
xseed3 = np.empty(0)
xseed4 = np.empty(0)

iis=0
jjs=0
kks=0
pps=0

for i in range(eseed):

    fnameref=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    optics_ref=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)
    optics_ref=optics_ref.reset_index(drop=True)

    #CORRECT: it 1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall1):
        xseed=np.append(xseed,iis)
        head_opt1=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt1)
        optics_all1=optics_all1.reset_index(drop=True)
        rms_disp_x=np.append(rms_disp_x,np.std((optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x=np.append(ave_disp_x,np.mean((optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y=np.append(rms_disp_y,np.std((optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y=np.append(ave_disp_y,np.mean((optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax3[0].plot(optics_all1['S']/1000., (optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]), ".")
        ax3[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax3[0].set_ylim(-0.05,0.05)
        ax3[1].plot(optics_all1['S']/1000., (optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]), ".")
        ax3[1].set_xlabel("longitundinal position [km]")
        ax3[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax3[1].set_ylim(-0.05,0.05)
        
        iis+=1
    else:
        continue
    
    #CORRECT: it 2
    fnameall2=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall2):
        xseed2=np.append(xseed2,jjs)
        head_opt2=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all2=pd.read_csv(fnameall2, skiprows=52, sep='\s+', names=head_opt2)
        optics_all2=optics_all2.reset_index(drop=True)
        rms_disp_x2=np.append(rms_disp_x2,np.std((optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x2=np.append(ave_disp_x2,np.mean((optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y2=np.append(rms_disp_y2,np.std((optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y2=np.append(ave_disp_y2,np.mean((optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax5[0].plot(optics_all2['S']/1000., (optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]), ".")
        ax5[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax5[0].set_ylim(-0.006,0.006)
        ax5[1].plot(optics_all2['S']/1000., (optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]), ".")
        ax5[1].set_xlabel("longitundinal position [km]")
        ax5[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax5[1].set_ylim(-0.002,0.002)
        
        jjs+=1
    else:
        continue


    #tune match: first it
    fnameall3=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall3):
        xseed3=np.append(xseed3,kks)
        head_opt3=pd.read_csv(fnameall3, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all3=pd.read_csv(fnameall3, skiprows=52, sep='\s+', names=head_opt3)
        optics_all3=optics_all3.reset_index(drop=True)
        rms_disp_x3=np.append(rms_disp_x3,np.std((optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x3=np.append(ave_disp_x3,np.mean((optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y3=np.append(rms_disp_y3,np.std((optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y3=np.append(ave_disp_y3,np.mean((optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax9[0].plot(optics_all3['S']/1000., (optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]), ".")
        ax9[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax9[0].set_ylim(-0.003,0.003)
        ax9[1].plot(optics_all3['S']/1000., (optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]), ".")
        ax9[1].set_xlabel("longitundinal position [km]")
        ax9[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax9[1].set_ylim(-0.003,0.003)
        
        kks+=1
    else:
        continue


    cnt=1
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(i+1,cnt)
    
    while os.path.exists(fnameall4):
        fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(i+1,cnt+1)
        cnt+=1
    
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs".format(i+1,cnt-1)

    #tune match: last it
    if os.path.exists(fnameall4):
        xseed4=np.append(xseed4,pps)
        head_opt4=pd.read_csv(fnameall4, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all4=pd.read_csv(fnameall4, skiprows=52, sep='\s+', names=head_opt4)
        optics_all4=optics_all4.reset_index(drop=True)
        rms_disp_x4=np.append(rms_disp_x4,np.std((optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x4=np.append(ave_disp_x4,np.mean((optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y4=np.append(rms_disp_y4,np.std((optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y4=np.append(ave_disp_y4,np.mean((optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax10[0].plot(optics_all4['S']/1000., (optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]), ".")
        ax10[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax10[0].set_ylim(-0.002,0.002) #for 90um
        #ax10[0].set_ylim(-0.003,0.003) #for 100um
        ax10[1].plot(optics_all4['S']/1000., (optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]), ".")
        ax10[1].set_xlabel("longitundinal position [km]")
        ax10[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax10[1].set_ylim(-0.002,0.002) #for 90um
        #ax10[1].set_ylim(-0.003,0.003) #for 100um
        
        pps+=1
    else:
        continue


#to heavy for pdf format
fig3.savefig(path+'disp_distribution_it1_{0}seeds.png'.format(eseed))
fig5.savefig(path+'disp_distribution_it2_{0}seeds.png'.format(eseed))
fig9.savefig(path+'disp_distribution_tune_match_first_it_{0}seeds.png'.format(eseed))
fig10.savefig(path+'disp_distribution_tune_match_last_it_{0}seeds.png'.format(eseed))

#plot of the rms: all
fig17, ax17=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax17[0].plot(xseed, rms_disp_x, ".", label='CORRECT it1')
ax17[0].plot(xseed2, rms_disp_x2, ".", label='CORRECT it2')
ax17[0].plot(xseed3, rms_disp_x3, ".", label='tune match first it')
ax17[0].plot(xseed4, rms_disp_x4, ".", label='tune match last it')
ax17[0].set_ylabel("rms$_x$ [m$^{1/2}$]")
ax17[0].legend(fontsize=8,loc='best')
ax17[1].plot(xseed, rms_disp_y, ".", label='CORRECT it1')
ax17[1].plot(xseed2, rms_disp_y2, ".", label='CORRECT it2')
ax17[1].plot(xseed3, rms_disp_y3, ".", label='tune match first it')
ax17[1].plot(xseed4, rms_disp_y4, ".", label='tune match last it')
ax17[1].set_xlabel("seed")
ax17[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax17[1].legend(fontsize=8,loc='best')
fig17, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig17.savefig(path+'rms_disp_all_{0}seeds.pdf'.format(eseed))

#plot of the mean: all
fig18, ax18=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax18[0].plot(xseed, ave_disp_x, ".", label='CORRECT it1')
ax18[0].plot(xseed2, ave_disp_x2, ".", label='CORRECT it2')
ax18[0].plot(xseed3, ave_disp_x3, ".", label='tune match first it')
ax18[0].plot(xseed4, ave_disp_x4, ".", label='tune match last it')
ax18[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax18[0].legend(fontsize=8,loc='best')
ax18[1].plot(xseed, ave_disp_y, ".", label='CORRECT it1')
ax18[1].plot(xseed2, ave_disp_y2, ".", label='CORRECT it2')
ax18[1].plot(xseed3, ave_disp_y3, ".", label='tune match first it')
ax18[1].plot(xseed4, ave_disp_y4, ".", label='tune match last it')
ax18[1].set_xlabel("seed")
ax18[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax18[1].legend(fontsize=8,loc='best')
fig18, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig18.savefig(path+'mean_disp_all_{0}seeds.pdf'.format(eseed))

#plot of the rms: it without tune match
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_disp_x, ".", label='it n°1')
ax2[0].plot(xseed2, rms_disp_x2, ".", label='it n°2')
ax2[0].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax2[0].legend(fontsize=10,loc='best')
ax2[1].plot(xseed, rms_disp_y, ".", label='it n°1')
ax2[1].plot(xseed2, rms_disp_y2, ".", label='it n°2')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax2[1].legend(fontsize=10,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_disp_correct_{0}seeds.pdf'.format(eseed))

#plot of the mean: it without tune match
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_disp_x, ".", label='it n°1')
ax4[0].plot(xseed2, ave_disp_x2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax4[0].legend(fontsize=10,loc='best')
ax4[1].plot(xseed, ave_disp_y, ".", label='it n°1')
ax4[1].plot(xseed2, ave_disp_y2, ".", label='it n°2')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax4[1].legend(fontsize=10,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_disp_correct_{0}seeds.pdf'.format(eseed))

#plot of the rms: it with tune match
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(xseed3, rms_disp_x3, ".", label='first it')
ax12[0].plot(xseed4, rms_disp_x4, ".", label='last it')
ax12[0].set_ylabel("rms$_x$ [m$^{1/2}$]")
ax12[0].legend(fontsize=10,loc='best')
ax12[1].plot(xseed3, rms_disp_y3, ".", label='first it')
ax12[1].plot(xseed4, rms_disp_y4, ".", label='last it')
ax12[1].set_xlabel("seed")
ax12[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax12[1].legend(fontsize=10,loc='best')
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig12.savefig(path+'rms_disp_tune_match_{0}seeds.pdf'.format(eseed))

#plot of the mean: it with tune match
fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax13[0].plot(xseed3, ave_disp_x3, ".", label='first it')
ax13[0].plot(xseed4, ave_disp_x4, ".", label='last it')
ax13[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax13[0].legend(fontsize=10,loc='best')
ax13[1].plot(xseed3, ave_disp_y3, ".", label='first it')
ax13[1].plot(xseed4, ave_disp_y4, ".", label='last it')
ax13[1].set_xlabel("seed")
ax13[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax13[1].legend(fontsize=10,loc='best')
fig13, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig13.savefig(path+'mean_disp_tune_match_{0}seeds.pdf'.format(eseed))

#for iteration n°1: it without tune match
fig1, ax1=plt.subplots(nrows=2, ncols=2, sharey=True)
ax1[0,0].hist(rms_disp_x,100)
ax1[0,0].set_xlabel("rms$_x$ [m$^{1/2}$]")
ax1[0,1].hist(ave_disp_x,100)
ax1[0,1].set_xlabel("mean$_x$ [m$^{1/2}$]")
ax1[0,0].set_ylabel("counts")
ax1[1,0].hist(rms_disp_y,100)
ax1[1,0].set_xlabel("rms$_y$ [m$^{1/2}$]")
ax1[1,1].hist(ave_disp_y,100)
ax1[1,1].set_xlabel("mean$_y$ [m$^{1/2}$]")
ax1[1,0].set_ylabel("counts")
fig1, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig1.tight_layout()
fig1.savefig(path+'hist_disp_correct_it1_{0}seeds.pdf'.format(eseed))

#for iteration n°2: it without tune match
fig8, ax8=plt.subplots(nrows=2, ncols=2, sharey=True)
ax8[0,0].hist(rms_disp_x2,100)
ax8[0,0].set_xlabel("rms$_x$ [m$^{1/2}$]")
ax8[0,1].hist(ave_disp_x2,100)
ax8[0,1].set_xlabel("mean$_x$ [m$^{1/2}$]")
ax8[0,0].set_ylabel("counts")
ax8[1,0].hist(rms_disp_y2,100)
ax8[1,0].set_xlabel("rms$_y$ [m$^{1/2}$]")
ax8[1,1].hist(ave_disp_y2,100)
ax8[1,1].set_xlabel("mean$_y$ [m$^{1/2}$]")
ax8[1,0].set_ylabel("counts")
fig8, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig8.tight_layout()
fig8.savefig(path+'hist_disp_correct_it2_{0}seeds.pdf'.format(eseed))

#for iteration n°1: it with tune match
fig14, ax14=plt.subplots(nrows=2, ncols=2, sharey=True)
ax14[0,0].hist(rms_disp_x3,100)
ax14[0,0].set_xlabel("rms$_x$ [m$^{1/2}$]")
ax14[0,1].hist(ave_disp_x3,100)
ax14[0,1].set_xlabel("mean$_x$ [m$^{1/2}$]")
ax14[0,0].set_ylabel("counts")
ax14[1,0].hist(rms_disp_y3,100)
ax14[1,0].set_xlabel("rms$_y$ [m$^{1/2}$]")
ax14[1,1].hist(ave_disp_y3,100)
ax14[1,1].set_xlabel("mean$_y$ [m$^{1/2}$]")
ax14[1,0].set_ylabel("counts")
fig14, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig14.tight_layout()
fig14.savefig(path+'hist_disp_tune_match_first_it_{0}seeds.pdf'.format(eseed))

#for iteration n°2: it with tune match
fig15, ax15=plt.subplots(nrows=2, ncols=2, sharey=True)
ax15[0,0].hist(rms_disp_x4,100)
ax15[0,0].set_xlabel("rms$_x$ [m$^{1/2}$]")
ax15[0,1].hist(ave_disp_x4,100)
ax15[0,1].set_xlabel("mean$_x$ [m$^{1/2}$]")
ax15[0,0].set_ylabel("counts")
ax15[1,0].hist(rms_disp_y4,100)
ax15[1,0].set_xlabel("rms$_y$ [m$^{1/2}$]")
ax15[1,1].hist(ave_disp_y4,100)
ax15[1,1].set_xlabel("mean$_y$ [m$^{1/2}$]")
ax15[1,0].set_ylabel("counts")
fig15, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig15.tight_layout()
fig15.savefig(path+'hist_disp_tune_match_last_it_{0}seeds.pdf'.format(eseed))




