#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the orbit after correction and tune match.
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
err_bpm=50 #bpm offset
eseed=100 #nb of seeds

#path definition
#path="/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_100seeds_corrhplus/".format(err_mq)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_dxdy_{0}_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_ms_dxdy_150_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_ms)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/'
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
fig21, ax21=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig21, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig22, ax22=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
fig22, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)

#line sextuoff it1
rms_dist_x=np.empty(0)
ave_dist_x=np.empty(0)
rms_dist_y=np.empty(0)
ave_dist_y=np.empty(0)

#line sextuoff it2
rms_dist_x2=np.empty(0)
ave_dist_x2=np.empty(0)
rms_dist_y2=np.empty(0)
ave_dist_y2=np.empty(0)

#ring sextuoff first it
rms_dist_x3=np.empty(0)
ave_dist_x3=np.empty(0)
rms_dist_y3=np.empty(0)
ave_dist_y3=np.empty(0)

#ring sextuoff last it
rms_dist_x4=np.empty(0)
ave_dist_x4=np.empty(0)
rms_dist_y4=np.empty(0)
ave_dist_y4=np.empty(0)

#ring sextuon
rms_dist_x5=np.empty(0)
ave_dist_x5=np.empty(0)
rms_dist_y5=np.empty(0)
ave_dist_y5=np.empty(0)

#tm successful it0
rms_dist_x6=np.empty(0)
ave_dist_x6=np.empty(0)
rms_dist_y6=np.empty(0)
ave_dist_y6=np.empty(0)

#tm all it0
rms_dist_x7 =np.empty(0)
ave_dist_x7 =np.empty(0)
rms_dist_y7 =np.empty(0)
ave_dist_y7 =np.empty(0)

#tm successful it1
rms_dist_x8 =np.empty(0)
ave_dist_x8 =np.empty(0)
rms_dist_y8 =np.empty(0)
ave_dist_y8 =np.empty(0)

#tm all it1
rms_dist_x9 =np.empty(0)
ave_dist_x9 =np.empty(0)
rms_dist_y9 =np.empty(0)
ave_dist_y9 =np.empty(0)

xseed = np.empty(0)
xseed2 = np.empty(0)
xseed3 = np.empty(0)
xseed4 = np.empty(0)
xseed5 = np.empty(0)
xseed6 = np.empty(0)
xseed7 = np.empty(0)
xseed8 = np.empty(0)
xseed9 = np.empty(0)

iis=0
jjs=0
kks=0
pps=0
qqs=0
tts=0
mms=0
ffs=0
dds=0

orbit_x_all_seed=np.empty(0)
orbit_y_all_seed=np.empty(0)

for j in range(eseed):
    
    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(j+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(j+1)
    mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)
    orbit_x_all_seed=np.append(orbit_x_all_seed,orbit_x)
    orbit_y_all_seed=np.append(orbit_y_all_seed,orbit_y)
    
max_orbit_x=np.max(orbit_x_all_seed)
max_orbit_y=np.max(orbit_y_all_seed)

print('analytical rms: max_orbit_x = ',max_orbit_x)
print('analytical rms: max_orbit_y = ',max_orbit_y)
print('\n')

for i in range(eseed):
    
    file1=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    file2=path+'test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(i+1)
    head_opt=pd.read_csv(file1, header=50, sep='\s+', nrows=0).columns[1:]
    headers=pd.read_csv(file1, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
    Q1=float(headers.iat[22,3]) #tune in plane x
    Q2=float(headers.iat[23,3]) #tune in plane y
    
    #line sextuoff it1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it1_seed{0}.tfs".format(i+1)
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
        if i==1:
            ax3[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax3[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax3[0].legend(fontsize=8,loc='best')
        #ax3[0].set_ylim(-250e-6,250e-6)
        ax3[0].set_ylim(-350e-6,350e-6)
        ax3[1].plot(optics_all1['S']/1000., optics_all1['Y'], ".")
        ax3[1].set_xlabel("longitundinal position [km]")
        ax3[1].set_ylabel("y [m]")
        if i==1:
            ax3[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax3[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax3[1].legend(fontsize=8,loc='best')
        #ax3[1].set_ylim(-250e-6,250e-6)
        ax3[1].set_ylim(-350e-6,350e-6)
        
        iis+=1
    else:
        continue
    
    #line sextuoff it2
    fnameall2=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs".format(i+1)
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
        if i==1:
            ax5[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax5[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax5[0].legend(fontsize=8,loc='best')
        #ax5[0].set_ylim(-250e-6,250e-6)
        ax5[0].set_ylim(-350e-6,350e-6)
        ax5[1].plot(optics_all2['S']/1000., optics_all2['Y'], ".")
        ax5[1].set_xlabel("longitundinal position [km]")
        ax5[1].set_ylabel("y [m]")
        if i==1:
            ax5[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax5[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax5[1].legend(fontsize=8,loc='best')
        #ax5[1].set_ylim(-250e-6,250e-6)
        ax5[1].set_ylim(-350e-6,350e-6)
        
        jjs+=1
    else:
        continue

    #ring sextuon
    fnameall5=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall5):
        xseed5=np.append(xseed5,qqs)
        head_opt5=pd.read_csv(fnameall5, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all5=pd.read_csv(fnameall5, skiprows=52, sep='\s+', names=head_opt5)
        optics_all5=optics_all5.reset_index(drop=True)
        rms_dist_x5=np.append(rms_dist_x5,np.std(optics_all5['X']))
        ave_dist_x5=np.append(ave_dist_x5,np.mean(optics_all5['X']))
        rms_dist_y5=np.append(rms_dist_y5,np.std(optics_all5['Y']))
        ave_dist_y5=np.append(ave_dist_y5,np.mean(optics_all5['Y']))

        ax11[0].plot(optics_all5['S']/1000., optics_all5['X'], ".")
        ax11[0].set_ylabel("x [m]")
        if i==1:
            ax11[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax11[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax11[0].legend(fontsize=8,loc='best')
        #ax11[0].set_ylim(-250e-6,250e-6)
        #ax11[0].set_ylim(-350e-6,350e-6) #200um
        #ax11[0].set_ylim(-500e-6,500e-6) #200um
        ax11[1].plot(optics_all5['S']/1000., optics_all5['Y'], ".")
        ax11[1].set_xlabel("longitundinal position [km]")
        ax11[1].set_ylabel("y [m]")
        if i==1:
            ax11[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax11[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax11[1].legend(fontsize=8,loc='best')
        #ax11[1].set_ylim(-250e-6,250e-6)
        #ax11[1].set_ylim(-350e-6,350e-6) #200um
        #ax11[1].set_ylim(-500e-6,500e-6) #200um
        
        qqs+=1
    else:
        continue

    #ring sextuoff first it
    fnameall3=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs".format(i+1)
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
        if i==1:
            ax9[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax9[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax9[0].legend(fontsize=8,loc='best')
        #ax9[0].set_ylim(-250e-6,250e-6)
        ax9[0].set_ylim(-350e-6,350e-6) #200um
        ax9[1].plot(optics_all3['S']/1000., optics_all3['Y'], ".")
        ax9[1].set_xlabel("longitundinal position [km]")
        ax9[1].set_ylabel("y [m]")
        if i==1:
            ax9[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax9[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax9[1].legend(fontsize=8,loc='best')
        #ax9[1].set_ylim(-250e-6,250e-6)
        ax9[1].set_ylim(-350e-6,350e-6) #200um       
        
        kks+=1
    else:
        continue

    cnt=1
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{1}_seed{0}.tfs".format(i+1,cnt)
    while os.path.exists(fnameall4):
        fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{1}_seed{0}.tfs".format(i+1,cnt+1)
        cnt+=1
        
    fnameall4=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{1}_seed{0}.tfs".format(i+1,cnt-1)
    
    #ring sextuoff last it
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
        if i==1:
            ax10[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax10[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax10[0].legend(fontsize=8,loc='best')
        #ax10[0].set_ylim(-250e-6,250e-6)
        ax10[0].set_ylim(-350e-6,350e-6) #200um
        ax10[1].plot(optics_all4['S']/1000., optics_all4['Y'], ".")
        ax10[1].set_xlabel("longitundinal position [km]")
        ax10[1].set_ylabel("y [m]")
        if i==1:
            ax10[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax10[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax10[1].legend(fontsize=8,loc='best')
        #ax10[1].set_ylim(-250e-6,250e-6)
        ax10[1].set_ylim(-350e-6,350e-6) #200um
        
        pps+=1
    else:
        continue
   
    #ring tm without svd
    fnameall6=path+"test_seed{0}/FCCee_heb_modett_tune_match_it0_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall6):
        head_opt=pd.read_csv(fnameall6, header=50, sep='\s+', nrows=0).columns[1:]
        headers=pd.read_csv(fnameall6, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
        Q_x=float(headers.iat[22,3]) #tune in plane x
        Q_y=float(headers.iat[23,3]) #tune in plane y

        head_opt6=pd.read_csv(fnameall6, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all6=pd.read_csv(fnameall6, skiprows=52, sep='\s+', names=head_opt6)
        optics_all6=optics_all6.reset_index(drop=True)

        if round(Q_x,3)==round(Q1,3) and round(Q_y,2)==round(Q2,2):
            xseed6=np.append(xseed6,tts)
            rms_dist_x6=np.append(rms_dist_x6,np.std(optics_all6['X']))
            ave_dist_x6=np.append(ave_dist_x6,np.mean(optics_all6['X']))
            rms_dist_y6=np.append(rms_dist_y6,np.std(optics_all6['Y']))
            ave_dist_y6=np.append(ave_dist_y6,np.mean(optics_all6['Y']))

            ax17[0].plot(optics_all6['S']/1000., optics_all6['X'], ".")
            ax17[0].set_ylabel("x [m]")
            if i==1:
                ax17[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
                ax17[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
                ax17[0].legend(fontsize=8,loc='best')
            #ax17[0].set_ylim(-250e-6,250e-6)
            #ax17[0].set_ylim(-350e-6,350e-6) #200um
            #ax17[0].set_ylim(-100e-5,100e-5) #200um
            ax17[1].plot(optics_all6['S']/1000., optics_all6['Y'], ".")
            ax17[1].set_xlabel("longitundinal position [km]")
            ax17[1].set_ylabel("y [m]")
            if i==1:
                ax17[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
                ax17[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
                ax17[1].legend(fontsize=8,loc='best')
            #ax17[1].set_ylim(-250e-6,250e-6)
            #ax17[1].set_ylim(-350e-6,350e-6) #200um
            #ax17[0].set_ylim(-100e-5,100e-5) #200um

            tts+=1

        xseed7=np.append(xseed7,mms)
        rms_dist_x7=np.append(rms_dist_x7,np.std(optics_all6['X']))
        ave_dist_x7=np.append(ave_dist_x7,np.mean(optics_all6['X']))
        rms_dist_y7=np.append(rms_dist_y7,np.std(optics_all6['Y']))
        ave_dist_y7=np.append(ave_dist_y7,np.mean(optics_all6['Y']))

        ax18[0].plot(optics_all6['S']/1000., optics_all6['X'], ".")
        ax18[0].set_ylabel("x [m]")
        if i==1:
            ax18[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax18[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax18[0].legend(fontsize=8,loc='best')
        #ax18[0].set_ylim(-250e-6,250e-6)
        #ax18[0].set_ylim(-350e-6,350e-6) #200um
        #ax18[0].set_ylim(-100e-5,100e-5) #200um
        ax18[1].plot(optics_all6['S']/1000., optics_all6['Y'], ".")
        ax18[1].set_xlabel("longitundinal position [km]")
        ax18[1].set_ylabel("y [m]")
        if i==1:
            ax18[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax18[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax18[1].legend(fontsize=8,loc='best')
        #ax18[1].set_ylim(-250e-6,250e-6)
        #ax18[1].set_ylim(-350e-6,350e-6) #200um
        #ax18[1].set_ylim(-100e-5,100e-5) #200um

        mms+=1
    else:
        continue

    #ring tm with svd
    fnameall7=path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall7):
        head_opt=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
        headers=pd.read_csv(fnameall7, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
        Q_x=float(headers.iat[22,3]) #tune in plane x
        Q_y=float(headers.iat[23,3]) #tune in plane y

        head_opt7=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all7=pd.read_csv(fnameall7, skiprows=52, sep='\s+', names=head_opt6)
        optics_all7=optics_all7.reset_index(drop=True)

        if round(Q_x,3)==round(Q1,3) and round(Q_y,2)==round(Q2,2):
            xseed8=np.append(xseed8,ffs)
            rms_dist_x8=np.append(rms_dist_x8,np.std(optics_all7['X']))
            ave_dist_x8=np.append(ave_dist_x8,np.mean(optics_all7['X']))
            rms_dist_y8=np.append(rms_dist_y8,np.std(optics_all7['Y']))
            ave_dist_y8=np.append(ave_dist_y8,np.mean(optics_all7['Y']))

            ax21[0].plot(optics_all7['S']/1000., optics_all7['X'], ".")
            ax21[0].set_ylabel("x [m]")
            if i==1:
                ax21[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
                ax21[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
                ax21[0].legend(fontsize=8,loc='best')
            #ax21[0].set_ylim(-250e-6,250e-6)
            #ax21[0].set_ylim(-350e-6,350e-6) #200um
            #ax21[0].set_ylim(-100e-5,100e-5) #200um
            ax21[1].plot(optics_all7['S']/1000., optics_all7['Y'], ".")
            ax21[1].set_xlabel("longitundinal position [km]")
            ax21[1].set_ylabel("y [m]")
            if i==1:
                ax21[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
                ax21[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
                ax21[1].legend(fontsize=8,loc='best')
            #ax21[1].set_ylim(-250e-6,250e-6)
            #ax21[1].set_ylim(-350e-6,350e-6) #200um
            #ax21[1].set_ylim(-100e-5,100e-5) #200um

            ffs+=1
    
        xseed8=np.append(xseed8,dds)
        rms_dist_x8=np.append(rms_dist_x8,np.std(optics_all7['X']))
        ave_dist_x8=np.append(ave_dist_x8,np.mean(optics_all7['X']))
        rms_dist_y8=np.append(rms_dist_y8,np.std(optics_all7['Y']))
        ave_dist_y8=np.append(ave_dist_y8,np.mean(optics_all7['Y']))

        ax22[0].plot(optics_all7['S']/1000., optics_all7['X'], ".")
        ax22[0].set_ylabel("x [m]")
        if i==1:
            ax22[0].axhline(y=3*max_orbit_x, color='r', linestyle='--', label='3 analytical rms')
            ax22[0].axhline(y=-3*max_orbit_x, color='r', linestyle='--')
            ax22[0].legend(fontsize=8,loc='best')
        #ax22[0].set_ylim(-250e-6,250e-6)
        #ax22[0].set_ylim(-350e-6,350e-6) #200um
        #ax22[0].set_ylim(-100e-5,100e-5) #200um
        ax22[1].plot(optics_all7['S']/1000., optics_all7['Y'], ".")
        ax22[1].set_xlabel("longitundinal position [km]")
        ax22[1].set_ylabel("y [m]")
        if i==1:
            ax22[1].axhline(y=3*max_orbit_y, color='r', linestyle='--', label='3 analytical rms')
            ax22[1].axhline(y=-3*max_orbit_y, color='r', linestyle='--')
            ax22[1].legend(fontsize=8,loc='best')
        #ax22[1].set_ylim(-250e-6,250e-6)
        #ax22[1].set_ylim(-350e-6,350e-6) #200um
        #ax22[1].set_ylim(-100e-5,100e-5) #200um

        dds+=1
    else:
        continue

print('rms_dist_x2[0] = ',rms_dist_x2[0])
print('rms_dist_y2[0] = ',rms_dist_y2[0])
print('\n')

mean_rms_dist_x2=np.mean(rms_dist_x2)
mean_rms_dist_y2=np.mean(rms_dist_y2)

print('numerical rms: mean_rms_dist_x2 = ',mean_rms_dist_x2)
print('numerical rms: mean_rms_dist_y2 = ',mean_rms_dist_y2)

#sys.exit()

#to heavy for pdf format
fig3.savefig(path+'orbit_distribution_line_sextuoff_it1_{0}seeds.png'.format(eseed+9))
fig5.savefig(path+'orbit_distribution_line_sextuoff_it2_{0}seeds.png'.format(eseed+9))
fig9.savefig(path+'orbit_distribution_ring_sextuoff_first_it_{0}seeds.png'.format(eseed+9))
fig10.savefig(path+'orbit_distribution_ring_sextuoff_last_it_{0}seeds.png'.format(eseed+9))
fig11.savefig(path+'orbit_distribution_ring_sextuon_it1_{0}seeds.png'.format(eseed+9))
fig17.savefig(path+'orbit_distribution_tm_seccessful_it0_{0}seeds.png'.format(eseed+9))
fig18.savefig(path+'orbit_distribution_tm_all_it0_{0}seeds.png'.format(eseed+9))
fig21.savefig(path+'orbit_distribution_tm_seccessful_it1_{0}seeds.png'.format(eseed+9))
fig22.savefig(path+'orbit_distribution_tm_all_it1_{0}seeds.png'.format(eseed+9))

#sys.exit()

#plot of the rms: sextuoff line
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_dist_x, ".", label='it n°1')
ax2[0].plot(xseed2, rms_dist_x2, ".", label='it n°2')
ax2[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax2[0].set_ylabel("rms$_x$ [m]")
#ax2[0].set_ylim(0,10e-5)
ax2[0].set_ylim(0,25e-5) #200um
ax2[0].legend(fontsize=8,loc='best')
ax2[1].plot(xseed, rms_dist_y, ".", label='it n°1')
ax2[1].plot(xseed2, rms_dist_y2, ".", label='it n°2')
ax2[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m]")
#ax2[1].set_ylim(0,10e-5)
ax2[1].set_ylim(0,25e-5) #200um
ax2[1].legend(fontsize=8,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_orbit_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the mean: sextuoff line
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_dist_x, ".", label='it n°1')
ax4[0].plot(xseed2, ave_dist_x2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [m]")
#ax4[0].set_ylim(-20e-7,20e-7)
ax4[0].set_ylim(-30e-7,30e-7) #200um
ax4[0].legend(fontsize=8,loc='best')
ax4[1].plot(xseed, ave_dist_y, ".", label='it n°1')
ax4[1].plot(xseed2, ave_dist_y2, ".", label='it n°2')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m]")
#ax4[1].set_ylim(-20e-7,20e-7)
ax4[1].set_ylim(-30e-7,30e-7) #200um
ax4[1].legend(fontsize=8,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_orbit_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the rms: sextuoff/on ring
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(xseed3, rms_dist_x3, ".", label='sextuoff: first it')
ax12[0].plot(xseed4, rms_dist_x4, ".", label='sextuoff: last it')
ax12[0].plot(xseed5, rms_dist_x5, ".", label='sextuon')
ax12[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax12[0].set_ylabel("rms$_x$ [m]")
#ax12[0].set_ylim(0,10e-5)
ax12[0].set_ylim(0,25e-5) #200um
ax12[0].legend(fontsize=8,loc='best')
ax12[1].plot(xseed3, rms_dist_y3, ".", label='sextuoff: first it')
ax12[1].plot(xseed4, rms_dist_y4, ".", label='sextuoff: last it')
ax12[1].plot(xseed5, rms_dist_y5, ".", label='sextuon')
ax12[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax12[1].set_xlabel("seed")
ax12[1].set_ylabel("rms$_y$ [m]")
#ax12[1].set_ylim(0,10e-5)
ax12[1].set_ylim(0,25e-5) #200um
ax12[1].legend(fontsize=8,loc='best')
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig12.savefig(path+'rms_orbit_ring_{0}seeds.pdf'.format(eseed))

#plot of the mean: sextuoff/on ring
fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax13[0].plot(xseed3, ave_dist_x3, ".", label='sextuoff: first it')
ax13[0].plot(xseed4, ave_dist_x4, ".", label='sextuoff: last it')
ax13[0].plot(xseed5, ave_dist_x5, ".", label='sextuon')
ax13[0].set_ylabel("mean$_x$ [m]")
#ax13[0].set_ylim(-20e-7,20e-7)
ax13[0].set_ylim(-30e-7,30e-7) #200um
ax13[0].legend(fontsize=8,loc='best')
ax13[1].plot(xseed3, ave_dist_y3, ".", label='sextuoff: first it')
ax13[1].plot(xseed4, ave_dist_y4, ".", label='sextuoff: last it')
ax13[1].plot(xseed5, ave_dist_x5, ".", label='sextuon')
ax13[1].set_xlabel("seed")
ax13[1].set_ylabel("mean$_y$ [m]")
#ax13[1].set_ylim(-20e-7,20e-7)
ax13[1].set_ylim(-30e-7,30e-7) #200um
ax13[1].legend(fontsize=8,loc='best')
fig13, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig13.savefig(path+'mean_orbit_ring_{0}seeds.pdf'.format(eseed))

#plot of the rms: tm without svd
fig20, ax20=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax20[0].plot(xseed6, rms_dist_x6, ".", label='tm success')
ax20[0].plot(xseed7, rms_dist_x7, ".", label='tm all')
ax20[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax20[0].set_ylabel("rms$_x$ [m]")
#ax20[0].set_ylim(0,10e-5)
ax20[0].set_ylim(0,25e-5) #200um
ax20[0].legend(fontsize=8,loc='best')
ax20[1].plot(xseed6, rms_dist_y6, ".", label='tm success')
ax20[1].plot(xseed7, rms_dist_y7, ".", label='tm all')
ax20[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax20[1].set_xlabel("seed")
ax20[1].set_ylabel("rms$_y$ [m]")
#ax20[1].set_ylim(0,10e-5)
ax20[1].set_ylim(0,25e-5) #200um
ax20[1].legend(fontsize=8,loc='best')
fig20, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig20.savefig(path+'rms_orbit_tm_it0_{0}seeds.pdf'.format(eseed))

#plot of the mean: tm without svd
fig19, ax19=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax19[0].plot(xseed6, ave_dist_x6, ".", label='tm success')
ax19[0].plot(xseed7, ave_dist_x7, ".", label='tm all')
ax19[0].set_ylabel("mean$_x$ [m]")
#ax19[0].set_ylim(-20e-7,20e-7)
ax19[0].set_ylim(-30e-7,30e-7) #200um
ax19[0].legend(fontsize=8,loc='best')
ax19[1].plot(xseed6, ave_dist_y6, ".", label='tm success')
ax19[1].plot(xseed7, ave_dist_y7, ".", label='tm all')
ax19[1].set_xlabel("seed")
ax19[1].set_ylabel("mean$_y$ [m]")
#ax19[1].set_ylim(-20e-7,20e-7)
ax19[1].set_ylim(-30e-7,30e-7) #200um
ax19[1].legend(fontsize=8,loc='best')
fig19, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig19.savefig(path+'mean_orbit_tm_it0_{0}seeds.pdf'.format(eseed))

#plot of the rms: tm with svd
fig23, ax23=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax23[0].plot(xseed6, rms_dist_x6, ".", label='tm success')
ax23[0].plot(xseed7, rms_dist_x7, ".", label='tm all')
ax23[0].axhline(y=max_orbit_x, color='r', linestyle='--', label='analytical rms')
ax23[0].set_ylabel("rms$_x$ [m]")
#ax23[0].set_ylim(0,10e-5)
ax23[0].set_ylim(0,25e-5) #200um
ax23[0].legend(fontsize=8,loc='best')
ax23[1].plot(xseed6, rms_dist_y6, ".", label='tm success')
ax23[1].plot(xseed7, rms_dist_y7, ".", label='tm all')
ax23[1].axhline(y=max_orbit_y, color='r', linestyle='--', label='analytical rms')
ax23[1].set_xlabel("seed")
ax23[1].set_ylabel("rms$_y$ [m]")
#ax23[1].set_ylim(0,10e-5)
ax23[1].set_ylim(0,25e-5) #200um
ax23[1].legend(fontsize=8,loc='best')
fig23, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig23.savefig(path+'rms_orbit_tm_it1_{0}seeds.pdf'.format(eseed))

#plot of the mean: tm with svd
fig24, ax24=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax24[0].plot(xseed6, ave_dist_x6, ".", label='tm success')
ax24[0].plot(xseed7, ave_dist_x7, ".", label='tm all')
ax24[0].set_ylabel("mean$_x$ [m]")
#ax24[0].set_ylim(-20e-7,20e-7)
ax24[0].set_ylim(-30e-7,30e-7) #200um
ax24[0].legend(fontsize=8,loc='best')
ax24[1].plot(xseed6, ave_dist_y6, ".", label='tm success')
ax24[1].plot(xseed7, ave_dist_y7, ".", label='tm all')
ax24[1].set_xlabel("seed")
ax24[1].set_ylabel("mean$_y$ [m]")
#ax24[1].set_ylim(-20e-7,20e-7)
ax24[1].set_ylim(-30e-7,30e-7) #200um
ax24[1].legend(fontsize=8,loc='best')
fig24, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig24.savefig(path+'mean_orbit_tm_it1_{0}seeds.pdf'.format(eseed))







'''

#line sextuoff it1
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
fig1.savefig(path+'hist_orbit_line_sextuoff_it1_{0}seeds.pdf'.format(eseed))

#line sextuoff it2
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
fig8.savefig(path+'hist_orbit_line_sextuoff_it2_{0}seeds.pdf'.format(eseed))

#ring sextuoff first it
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
fig14.savefig(path+'hist_orbit_ring_sextuoff_first_it_{0}seeds.pdf'.format(eseed))

#ring sextuoff last it
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
fig15.savefig(path+'hist_orbit_ring_sextuoff_last_it_{0}seeds.pdf'.format(eseed))

#ring sextuon
fig16, ax16=plt.subplots(nrows=2, ncols=2, sharey=True)
ax16[0,0].hist(rms_dist_x5,100)
ax16[0,0].set_xlabel("rms$_x$ [m]")
ax16[0,1].hist(ave_dist_x5,100)
ax16[0,1].set_xlabel("mean$_x$ [m]")
ax16[0,0].set_ylabel("counts")
ax16[1,0].hist(rms_dist_y5,100)
ax16[1,0].set_xlabel("rms$_y$ [m]")
ax16[1,1].hist(ave_dist_y5,100)
ax16[1,1].set_xlabel("mean$_y$ [m]")
ax16[1,0].set_ylabel("counts")
fig16, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig16.tight_layout()
fig16.savefig(path+'hist_orbit_ring_sextuon_it1_{0}seeds.pdf'.format(eseed))

'''
