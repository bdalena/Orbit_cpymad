#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

'''
This python file calculates and saves the statistical study of the normalized dispersion after correction  and tune match.
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

err_mq=60 #quads offset
err_ms=200 #sextupoles offset
err_bpm=200 #bpm offset
eseed=100 #nb of seeds

#path definition
#path="/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_100seeds_corrhplus/".format(err_mq)
path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_dxdy_{0}_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_ms)
#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_ms_dxdy_150_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/'.format(err_bpm)

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
rms_disp_x =np.empty(0)
ave_disp_x =np.empty(0)
rms_disp_y =np.empty(0)
ave_disp_y =np.empty(0)

#line sextuoff it2
rms_disp_x2 =np.empty(0)
ave_disp_x2 =np.empty(0)
rms_disp_y2 =np.empty(0)
ave_disp_y2 =np.empty(0)

#ring sextuoff first it
rms_disp_x3 =np.empty(0)
ave_disp_x3 =np.empty(0)
rms_disp_y3 =np.empty(0)
ave_disp_y3 =np.empty(0)

#ring sextuoff last it
rms_disp_x4 =np.empty(0)
ave_disp_x4 =np.empty(0)
rms_disp_y4 =np.empty(0)
ave_disp_y4 =np.empty(0)

#ring sextuon
rms_disp_x5 =np.empty(0)
ave_disp_x5 =np.empty(0)
rms_disp_y5 =np.empty(0)
ave_disp_y5 =np.empty(0)

#tm successful it0
rms_disp_x6=np.empty(0)
ave_disp_x6=np.empty(0)
rms_disp_y6=np.empty(0)
ave_disp_y6=np.empty(0)

#tm all it0
rms_disp_x7 =np.empty(0)
ave_disp_x7 =np.empty(0)
rms_disp_y7 =np.empty(0)
ave_disp_y7 =np.empty(0)

#tm successful it1
rms_disp_x8 =np.empty(0)
ave_disp_x8 =np.empty(0)
rms_disp_y8 =np.empty(0)
ave_disp_y8 =np.empty(0)

#tm all it1
rms_disp_x9 =np.empty(0)
ave_disp_x9 =np.empty(0)
rms_disp_y9 =np.empty(0)
ave_disp_y9 =np.empty(0)

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

for i in range(eseed):

    fnameref=path+'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(i+1)
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    headers=pd.read_csv(fnameref, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
    Q1=float(headers.iat[22,3]) #tune in plane x
    Q2=float(headers.iat[23,3]) #tune in plane y
    optics_ref=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)
    optics_ref=optics_ref.reset_index(drop=True)

    #line sextuoff it1
    fnameall1=path+"test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall1):
        xseed=np.append(xseed,iis)
        head_opt1=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt1)
        optics_all1=optics_all1.reset_index(drop=True)
        rms_disp_x=np.append(rms_disp_x,np.std((optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x=np.append(ave_disp_x,np.mean((optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y=np.append(rms_disp_y,np.std((optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y=np.append(ave_disp_y,np.mean((optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax3[0].plot(optics_all1['S']/1000., (optics_all1["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax3[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax3[0].set_ylim(-0.01,0.01)
        ax3[1].plot(optics_all1['S']/1000., (optics_all1["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax3[1].set_xlabel("longitundinal position [km]")
        ax3[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax3[1].set_ylim(-0.01,0.01)
        
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
        rms_disp_x2=np.append(rms_disp_x2,np.std((optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x2=np.append(ave_disp_x2,np.mean((optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y2=np.append(rms_disp_y2,np.std((optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y2=np.append(ave_disp_y2,np.mean((optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax5[0].plot(optics_all2['S']/1000., (optics_all2["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax5[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax5[0].set_ylim(-0.002,0.002)
        ax5[1].plot(optics_all2['S']/1000., (optics_all2["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax5[1].set_xlabel("longitundinal position [km]")
        ax5[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax5[1].set_ylim(-0.002,0.002)
        
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
        rms_disp_x5=np.append(rms_disp_x5,np.std((optics_all5["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x5=np.append(ave_disp_x5,np.mean((optics_all5["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y5=np.append(rms_disp_y5,np.std((optics_all5["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y5=np.append(ave_disp_y5,np.mean((optics_all5["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax11[0].plot(optics_all5['S']/1000., (optics_all5["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax11[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        #ax11[0].set_ylim(-0.002,0.002)
        #ax11[0].set_ylim(-0.004,0.004) #200um
        ax11[1].plot(optics_all5['S']/1000., (optics_all5["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax11[1].set_xlabel("longitundinal position [km]")
        ax11[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        #ax11[1].set_ylim(-0.002,0.002)
        #ax11[1].set_ylim(-0.004,0.004) #200um
        
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
        rms_disp_x3=np.append(rms_disp_x3,np.std((optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x3=np.append(ave_disp_x3,np.mean((optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y3=np.append(rms_disp_y3,np.std((optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y3=np.append(ave_disp_y3,np.mean((optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax9[0].plot(optics_all3['S']/1000., (optics_all3["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax9[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        ax9[0].set_ylim(-0.06,0.06)
        ax9[1].plot(optics_all3['S']/1000., (optics_all3["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax9[1].set_xlabel("longitundinal position [km]")
        ax9[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax9[1].set_ylim(-0.004,0.004)
        
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
        rms_disp_x4=np.append(rms_disp_x4,np.std((optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x4=np.append(ave_disp_x4,np.mean((optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y4=np.append(rms_disp_y4,np.std((optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y4=np.append(ave_disp_y4,np.mean((optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax10[0].plot(optics_all4['S']/1000., (optics_all4["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax10[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        #ax10[0].set_ylim(-0.02,0.02)
        ax10[0].set_ylim(-0.04,0.04) #200um
        ax10[1].plot(optics_all4['S']/1000., (optics_all4["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax10[1].set_xlabel("longitundinal position [km]")
        ax10[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        ax10[1].set_ylim(-0.004,0.004)
        
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
            head_opt6=pd.read_csv(fnameall6, header=50, sep='\s+', nrows=0).columns[1:]
            optics_all6=pd.read_csv(fnameall6, skiprows=52, sep='\s+', names=head_opt6)
            optics_all6=optics_all6.reset_index(drop=True)
            rms_disp_x6=np.append(rms_disp_x6,np.std((optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
            ave_disp_x6=np.append(ave_disp_x6,np.mean((optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
            rms_disp_y6=np.append(rms_disp_y6,np.std((optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
            ave_disp_y6=np.append(ave_disp_y6,np.mean((optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

            ax17[0].plot(optics_all6['S']/1000., (optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
            ax17[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
            #ax17[0].set_ylim(-0.02,0.02)
            #ax17[0].set_ylim(-0.004,0.004) #200um
            ax17[1].plot(optics_all6['S']/1000., (optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
            ax17[1].set_xlabel("longitundinal position [km]")
            ax17[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
            #ax17[1].set_ylim(-0.004,0.004)

            tts+=1

        xseed7=np.append(xseed7,mms)
        rms_disp_x7=np.append(rms_disp_x7,np.std((optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x7=np.append(ave_disp_x7,np.mean((optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y7=np.append(rms_disp_y7,np.std((optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y7=np.append(ave_disp_y7,np.mean((optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax18[0].plot(optics_all6['S']/1000., (optics_all6["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax18[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        #ax18[0].set_ylim(-0.02,0.02)
        #ax18[0].set_ylim(-0.004,0.004) #200um
        ax18[1].plot(optics_all6['S']/1000., (optics_all6["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax18[1].set_xlabel("longitundinal position [km]")
        ax18[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        #ax18[1].set_ylim(-0.004,0.004)

        mms+=1


    #ring tm with svd
    fnameall7=path+"test_seed{0}/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(i+1)
    if os.path.exists(fnameall7):
        head_opt=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
        headers=pd.read_csv(fnameall7, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
        Q_x=float(headers.iat[22,3]) #tune in plane x
        Q_y=float(headers.iat[23,3]) #tune in plane y

        head_opt7=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
        optics_all7=pd.read_csv(fnameall7, skiprows=52, sep='\s+', names=head_opt7)
        optics_all7=optics_all7.reset_index(drop=True)

        if round(Q_x,3)==round(Q1,3) and round(Q_y,2)==round(Q2,2):
            xseed8=np.append(xseed8,ffs)
            head_opt7=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
            optics_all7=pd.read_csv(fnameall7, skiprows=52, sep='\s+', names=head_opt7)
            optics_all7=optics_all7.reset_index(drop=True)
            rms_disp_x8=np.append(rms_disp_x8,np.std((optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
            ave_disp_x8=np.append(ave_disp_x8,np.mean((optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
            rms_disp_y8=np.append(rms_disp_y8,np.std((optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
            ave_disp_y8=np.append(ave_disp_y8,np.mean((optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

            ax21[0].plot(optics_all7['S']/1000., (optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
            ax21[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
            #ax21[0].set_ylim(-0.02,0.02)
            #ax21[0].set_ylim(-0.004,0.004) #200um
            ax21[1].plot(optics_all7['S']/1000., (optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
            ax21[1].set_xlabel("longitundinal position [km]")
            ax21[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
            #ax21[1].set_ylim(-0.004,0.004)

            ffs+=1

        xseed9=np.append(xseed9,dds)
        rms_disp_x9=np.append(rms_disp_x9,np.std((optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        ave_disp_x9=np.append(ave_disp_x9,np.mean((optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"])))
        rms_disp_y9=np.append(rms_disp_y9,np.std((optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))
        ave_disp_y9=np.append(ave_disp_y9,np.mean((optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"])))

        ax22[0].plot(optics_all7['S']/1000., (optics_all7["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]))
        ax22[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
        #ax22[0].set_ylim(-0.02,0.02)
        #ax22[0].set_ylim(-0.004,0.004) #200um
        ax22[1].plot(optics_all7['S']/1000., (optics_all7["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]))
        ax22[1].set_xlabel("longitundinal position [km]")
        ax22[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
        #ax22[1].set_ylim(-0.004,0.004)

        dds+=1



#to heavy for pdf format
fig3.savefig(path+'disp_distribution_line_sextuoff_it1_{0}seeds.png'.format(eseed))
fig5.savefig(path+'disp_distribution_line_sextuoff_it2_{0}seeds.png'.format(eseed))
fig9.savefig(path+'disp_distribution_ring_sextuoff_first_it_{0}seeds.png'.format(eseed))
fig10.savefig(path+'disp_distribution_ring_sextuoff_last_it_{0}seeds.png'.format(eseed))
fig11.savefig(path+'disp_distribution_ring_sextuon_it1_{0}seeds.png'.format(eseed))
fig17.savefig(path+'disp_distribution_tm_seccessful_it0_{0}seeds.png'.format(eseed))
fig18.savefig(path+'disp_distribution_tm_all_it0_{0}seeds.png'.format(eseed))
fig21.savefig(path+'disp_distribution_tm_seccessful_it1_{0}seeds.png'.format(eseed))
fig22.savefig(path+'disp_distribution_tm_all_it1_{0}seeds.png'.format(eseed))

#plot of the rms: all
fig17, ax17=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax17[0].plot(xseed, rms_disp_x, ".", label='MS_off_line: it1')
ax17[0].plot(xseed2, rms_disp_x2, ".", label='MS_off_line: it2')
ax17[0].plot(xseed3, rms_disp_x3, ".", label='MS_off_ring: first it')
ax17[0].plot(xseed4, rms_disp_x4, ".", label='MS_off_ring: last it')
ax17[0].plot(xseed5, rms_disp_x5, ".", label='MS_on')
ax17[0].set_ylabel("rms$_x$ [m$^{1/2}$]")
ax17[0].legend(fontsize=8,loc='best')
ax17[1].plot(xseed, rms_disp_y, ".", label='MS_off_line: it1')
ax17[1].plot(xseed2, rms_disp_y2, ".", label='MS_off_line: it2')
ax17[1].plot(xseed3, rms_disp_y3, ".", label='MS_off_ring: first it')
ax17[1].plot(xseed4, rms_disp_y4, ".", label='MS_off_ring: last it')
ax17[1].plot(xseed5, rms_disp_y5, ".", label='MS_on')
ax17[1].set_xlabel("seed")
ax17[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax17[1].legend(fontsize=8,loc='best')
fig17, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig17.savefig(path+'rms_disp_all_{0}seeds.pdf'.format(eseed))

#plot of the mean: all
fig18, ax18=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax18[0].plot(xseed, ave_disp_x, ".", label='MS_off_line: it1')
ax18[0].plot(xseed2, ave_disp_x2, ".", label='MS_off_line: it2')
ax18[0].plot(xseed3, ave_disp_x3, ".", label='MS_off_ring: first it')
ax18[0].plot(xseed4, ave_disp_x4, ".", label='MS_off_ring: last it')
ax18[0].plot(xseed5, ave_disp_x5, ".", label='MS_on')
ax18[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax18[0].legend(fontsize=8,loc='best')
ax18[1].plot(xseed, ave_disp_y, ".", label='MS_off_line: it1')
ax18[1].plot(xseed2, ave_disp_y2, ".", label='MS_off_line: it2')
ax18[1].plot(xseed3, ave_disp_y3, ".", label='MS_off_ring: first it')
ax18[1].plot(xseed4, ave_disp_y4, ".", label='MS_off_ring: last it')
ax18[1].plot(xseed5, ave_disp_y5, ".", label='MS_on')
ax18[1].set_xlabel("seed")
ax18[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax18[1].legend(fontsize=8,loc='best')
fig18, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig18.savefig(path+'mean_disp_all_{0}seeds.pdf'.format(eseed))

#plot of the rms: sextuoff line
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(xseed, rms_disp_x, ".", label='it n°1')
ax2[0].plot(xseed2, rms_disp_x2, ".", label='it n°2')
ax2[0].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax2[0].legend(fontsize=8,loc='best')
ax2[1].plot(xseed, rms_disp_y, ".", label='it n°1')
ax2[1].plot(xseed2, rms_disp_y2, ".", label='it n°2')
ax2[1].set_xlabel("seed")
ax2[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax2[1].legend(fontsize=8,loc='best')
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+'rms_disp_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the mean: sextuoff line
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(xseed, ave_disp_x, ".", label='it n°1')
ax4[0].plot(xseed2, ave_disp_x2, ".", label='it n°2')
ax4[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax4[0].legend(fontsize=8,loc='best')
ax4[1].plot(xseed, ave_disp_y, ".", label='it n°1')
ax4[1].plot(xseed2, ave_disp_y2, ".", label='it n°2')
ax4[1].set_xlabel("seed")
ax4[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax4[1].legend(fontsize=8,loc='best')
fig4, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig4.savefig(path+'mean_disp_line_sextuoff_{0}seeds.pdf'.format(eseed))

#plot of the rms: sextuoff/on ring
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(xseed3, rms_disp_x3, ".", label='MS_off: first it')
ax12[0].plot(xseed4, rms_disp_x4, ".", label='MS_off: last it')
ax12[0].plot(xseed5, rms_disp_x5, ".", label='MS_on')
ax12[0].set_ylabel("rms$_x$ [m$^{1/2}$]")
ax12[0].legend(fontsize=8,loc='best')
ax12[1].plot(xseed3, rms_disp_y3, ".", label='MS_off: first it')
ax12[1].plot(xseed4, rms_disp_y4, ".", label='MS_off: last it')
ax12[1].plot(xseed5, rms_disp_y5, ".", label='MS_on')
ax12[1].set_xlabel("seed")
ax12[1].set_ylabel("rms$_y$ [m$^{1/2}$]")
ax12[1].legend(fontsize=8,loc='best')
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig12.savefig(path+'rms_disp_ring_{0}seeds.pdf'.format(eseed))

#plot of the mean: sextuoff/on ring
fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax13[0].plot(xseed3, ave_disp_x3, ".", label='MS_off: first it')
ax13[0].plot(xseed4, ave_disp_x4, ".", label='MS_off: last it')
ax13[0].plot(xseed4, ave_disp_x4, ".", label='MS_on')
ax13[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax13[0].legend(fontsize=8,loc='best')
ax13[1].plot(xseed3, ave_disp_y3, ".", label='MS_off: first it')
ax13[1].plot(xseed4, ave_disp_y4, ".", label='MS_off: last it')
ax13[1].plot(xseed4, ave_disp_y4, ".", label='MS_on')
ax13[1].set_xlabel("seed")
ax13[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax13[1].legend(fontsize=8,loc='best')
fig13, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig13.savefig(path+'mean_disp_ring_{0}seeds.pdf'.format(eseed))

#plot of the rms: tm without svd
fig20, ax20=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax20[0].plot(xseed6, rms_disp_x6, ".", label='tm success')
ax20[0].plot(xseed7, rms_disp_x7, ".", label='tm all')
ax20[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax20[0].legend(fontsize=8,loc='best')
ax20[1].plot(xseed6, rms_disp_y6, ".", label='tm success')
ax20[1].plot(xseed7, rms_disp_y7, ".", label='tm all')
ax20[1].set_xlabel("seed")
ax20[1].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax20[1].legend(fontsize=8,loc='best')
fig20, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig20.savefig(path+'rms_disp_tm_it0_{0}seeds.pdf'.format(eseed))

#plot of the mean: tm without svd
fig19, ax19=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax19[0].plot(xseed6, ave_disp_x6, ".", label='tm success')
ax19[0].plot(xseed7, ave_disp_x7, ".", label='tm all')
ax19[0].set_ylabel("mean$_x$ [m$^{1/2}$]")
ax19[0].legend(fontsize=8,loc='best')
ax19[1].plot(xseed6, ave_disp_y6, ".", label='tm success')
ax19[1].plot(xseed7, ave_disp_y7, ".", label='tm all')
ax19[1].set_xlabel("seed")
ax19[1].set_ylabel("mean$_y$ [m$^{1/2}$]")
ax19[1].legend(fontsize=8,loc='best')
fig19, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig19.savefig(path+'mean_disp_tm_it0_{0}seeds.pdf'.format(eseed))

#plot of the rms: tm with svd
fig23, ax23=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax23[0].plot(xseed8, rms_disp_x8, ".", label='tm success')
ax23[0].plot(xseed9, rms_disp_x9, ".", label='tm all')
ax23[0].set_ylabel("mean$_x$ [%]")
ax23[0].legend(fontsize=8,loc='best')
ax23[1].plot(xseed8, rms_disp_y8, ".", label='tm success')
ax23[1].plot(xseed9, rms_disp_y9, ".", label='tm all')
ax23[1].set_xlabel("seed")
ax23[1].set_ylabel("mean$_x$ [%]")
ax23[1].legend(fontsize=8,loc='best')
fig23, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig23.savefig(path+'rms_disp_tm_it1_{0}seeds.pdf'.format(eseed))

#plot of the mean: tm with svd
fig24, ax24=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax24[0].plot(xseed6, ave_disp_x6, ".", label='tm success')
ax24[0].plot(xseed7, ave_disp_x7, ".", label='tm all')
ax24[0].set_ylabel("mean$_x$ [%]")
ax24[0].legend(fontsize=8,loc='best')
ax24[1].plot(xseed6, ave_disp_y6, ".", label='tm success')
ax24[1].plot(xseed7, ave_disp_y7, ".", label='tm all')
ax24[1].set_xlabel("seed")
ax24[1].set_ylabel("mean$_y$ [%]")
ax24[1].legend(fontsize=8,loc='best')
fig24, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig24.savefig(path+'mean_disp_tm_it1_{0}seeds.pdf'.format(eseed))




'''
#line sextuoff it1
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

#line sextuoff it2
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

#ring sextuoff first it
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

#ring sextuoff last it
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


#ring sextuon
fig16, ax16=plt.subplots(nrows=2, ncols=2, sharey=True)
ax16[0,0].hist(rms_disp_x5,100)
ax16[0,0].set_xlabel("rms$_x$ [m$^{1/2}$]")
ax16[0,1].hist(ave_disp_x5,100)
ax16[0,1].set_xlabel("mean$_x$ [m$^{1/2}$]")
ax16[0,0].set_ylabel("counts")
ax16[1,0].hist(rms_disp_y5,100)
ax16[1,0].set_xlabel("rms$_y$ [m$^{1/2}$]")
ax16[1,1].hist(ave_disp_y5,100)
ax16[1,1].set_xlabel("mean$_y$ [m$^{1/2}$]")
ax16[1,0].set_ylabel("counts")
fig16, plt.subplots_adjust(left=.12, right=.97, top=.97, bottom=.10)
fig16.tight_layout()
fig16.savefig(path+'hist_orbit_ring_sextuon_it1_{0}seeds.pdf'.format(eseed))
'''
