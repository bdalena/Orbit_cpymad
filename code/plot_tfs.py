#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

#-----DATAFRAME-----

plt.rcParams.update({'font.size':13})
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

err_mq=150
err_ms=150
err_bpm=150
res_bpm=50
eseed=3

#path="/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_100seeds_corrhplus/test_seed{1}".format(err_mq,eseed)

#path='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_dxdy_{0}_mb_b1r_dpsi_mq_dxdy_150_tm_100seeds/test_seed{1}'.format(err_bpm,eseed)

#path ='/feynman/work/dacm/leda/td271008/Orbit_cpymad/mb_b1r_dpsi_mq_dxdy_{0}_tm_100seeds/test_seed{1}'.format(err_mq,eseed)

path ='/feynman/work/dacm/leda/td271008/Orbit_cpymad/ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/test_seed{1}'.format(err_ms,eseed)

#path ='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_res_{0}_ms_dxdy_150_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/test_seed{1}'.format(res_bpm,eseed)

#path ='/feynman/work/dacm/leda/td271008/Orbit_cpymad/bpm_dxdy_{0}_mb_b1r_dpsi_mq_dxdy_150_tm_svd_100seeds/test_seed{1}'.format(err_bpm,eseed)


fnameref=path+"/FCCee_heb_modett_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
optics_ref=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

fnamerr=path+"/FCCee_heb_modett_err_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnamerr, header=50, sep='\s+', nrows=0).columns[1:]
optics_err=pd.read_csv(fnamerr, skiprows=52, sep='\s+', names=head_opt)
optics_err=optics_err.reset_index(drop=True)

fnametre=path+"/FCCee_heb_modett_orbcor_arc4.tfs".format(eseed) # NB: if cycling on an other IP than IP5, change the file for the last one done by the teading
head_opt=pd.read_csv(fnametre, header=50, sep='\s+', nrows=0).columns[1:]
optics_tre=pd.read_csv(fnametre, skiprows=52, sep='\s+', names=head_opt)
optics_tre=optics_tre.reset_index(drop=True)

#-----PLOT-----

# Plot orbits in m
fig1, ax1=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax1[0].plot(optics_err['S']/1000., optics_err['X'], "-r")
ax1[0].set_ylabel("x [m]")
ax1[1].plot(optics_err['S']/1000., optics_err['Y'], "-b")
ax1[1].set_xlabel("longitundinal position [km]")
ax1[1].set_ylabel("y [m]")
fig1, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig1.savefig(path+"/orbit_line_err_seed{0}.pdf".format(eseed))

fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(optics_tre['S']/1000., optics_tre['X'], "-r")
ax2[0].set_ylabel("x [m]")
ax2[1].plot(optics_tre['S']/1000., optics_tre['Y'], "-b")
ax2[1].set_xlabel("longitundinal position [km]")
ax2[1].set_ylabel("y [m]")
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+"/orbit_line_err_treading_seed{0}.pdf".format(eseed))


if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_line_sextuoff_it1_seed{0}.tfs".format(eseed)):
    fnameall11=path+"/FCCee_heb_modett_orbcor_all_line_sextuoff_it1_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall11, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all11=pd.read_csv(fnameall11, skiprows=52, sep='\s+', names=head_opt)
    optics_all11=optics_all11.reset_index(drop=True)
    fig31, ax31=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax31[0].plot(optics_all11['S']/1000., optics_all11['X'], "-r")
    ax31[0].set_ylabel("x [m]")
    ax31[1].plot(optics_all11['S']/1000., optics_all11['Y'], "-b")
    ax31[1].set_xlabel("longitundinal position [km]")
    ax31[1].set_ylabel("y [m]")
    fig31, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig31.savefig(path+"/orbit_line_err_treading_correctall_sextuoff_it1_seed{0}.pdf".format(eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs".format(eseed)):
    fnameall21=path+"/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs".format(eseed) 
    head_opt=pd.read_csv(fnameall21, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all21=pd.read_csv(fnameall21, skiprows=52, sep='\s+', names=head_opt)
    optics_all21=optics_all21.reset_index(drop=True)
    fig41, ax41=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax41[0].plot(optics_all21['S']/1000., optics_all21['X'], "-r")
    ax41[0].set_ylabel("x [m]")
    ax41[1].plot(optics_all21['S']/1000., optics_all21['Y'], "-b")
    ax41[1].set_xlabel("longitundinal position [km]")
    ax41[1].set_ylabel("y [m]")
    fig41, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig41.savefig(path+"/orbit_line_err_treading_correctall_sextuoff_it2_seed{0}.pdf".format(eseed))


if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs".format(eseed)):
    fnameall22=path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs".format(eseed) 
    head_opt=pd.read_csv(fnameall22, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all22=pd.read_csv(fnameall22, skiprows=52, sep='\s+', names=head_opt)
    optics_all22=optics_all22.reset_index(drop=True)
    fig42, ax42=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax42[0].plot(optics_all22['S']/1000., optics_all22['X'], "-r")
    ax42[0].set_ylabel("x [m]")
    ax42[1].plot(optics_all22['S']/1000., optics_all22['Y'], "-b")
    ax42[1].set_xlabel("longitundinal position [km]")
    ax42[1].set_ylabel("y [m]")
    fig42, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig42.savefig(path+"/orbit_ring_err_treading_correctall_sextuoff_first_it_seed{0}.pdf".format(eseed))

    
'''
cnt=1
fnameall23=path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{0}_seed{1}.tfs".format(cnt,eseed)
while os.path.exists(fnameall23):
    fnameall23=path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{0}_seed{1}.tfs".format(cnt+1,eseed)
    cnt+=1

fnameall23=path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{0}_seed{1}.tfs".format(cnt-1,eseed)
print(cnt)
print(fnameall23)
'''

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it8_seed{0}.tfs".format(eseed)):
    fnameall23=path+"/FCCee_heb_modett_orbcor_all_ring_sextuoff_it8_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall23, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all23=pd.read_csv(fnameall23, skiprows=52, sep='\s+', names=head_opt)
    optics_all23=optics_all23.reset_index(drop=True)
    fig43, ax43=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax43[0].plot(optics_all23['S']/1000., optics_all23['X'], "-r")
    ax43[0].set_ylabel("x [m]")
    ax43[1].plot(optics_all23['S']/1000., optics_all23['Y'], "-b")
    ax43[1].set_xlabel("longitundinal position [km]")
    ax43[1].set_ylabel("y [m]")
    fig43, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig43.savefig(path+"/orbit_ring_err_treading_correctall_sextuoff_last_it_seed{0}.pdf".format(eseed))


if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)):
    fnameall12=path+"/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall12, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all12=pd.read_csv(fnameall12, skiprows=52, sep='\s+', names=head_opt)
    optics_all12=optics_all12.reset_index(drop=True)
    fig32, ax32=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax32[0].plot(optics_all12['S']/1000., optics_all12['X'], "-r")
    ax32[0].set_ylabel("x [m]")
    ax32[1].plot(optics_all12['S']/1000., optics_all12['Y'], "-b")
    ax32[1].set_xlabel("longitundinal position [km]")
    ax32[1].set_ylabel("y [m]")
    fig32, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig32.savefig(path+"/orbit_ring_err_treading_correctall_sextuon_it1_seed{0}.pdf".format(eseed))


if os.path.exists(path+"/FCCee_heb_modett_tune_match_it0_seed{0}.tfs".format(eseed)):
    fnameall13=path+"/FCCee_heb_modett_tune_match_it0_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall13, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all13=pd.read_csv(fnameall13, skiprows=52, sep='\s+', names=head_opt)
    optics_all13=optics_all13.reset_index(drop=True)
    fig33, ax33=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax33[0].plot(optics_all13['S']/1000., optics_all13['X'], "-r")
    ax33[0].set_ylabel("x [m]")
    ax33[1].plot(optics_all13['S']/1000., optics_all13['Y'], "-b")
    ax33[1].set_xlabel("longitundinal position [km]")
    ax33[1].set_ylabel("y [m]")
    fig33, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig33.savefig(path+"/orbit_ring_err_treading_correctall_tm_it0_seed{0}.pdf".format(eseed))

if os.path.exists(path+"/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(eseed)):
    fnameall14=path+"/FCCee_heb_modett_tune_match_it1_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall14, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all14=pd.read_csv(fnameall14, skiprows=52, sep='\s+', names=head_opt)
    optics_all14=optics_all14.reset_index(drop=True)
    fig34, ax34=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax34[0].plot(optics_all14['S']/1000., optics_all14['X'], "-r")
    ax34[0].set_ylabel("x [m]")
    ax34[1].plot(optics_all14['S']/1000., optics_all14['Y'], "-b")
    ax34[1].set_xlabel("longitundinal position [km]")
    ax34[1].set_ylabel("y [m]")
    fig34, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig34.savefig(path+"/orbit_ring_err_treading_correctall_tm_it1_seed{0}.pdf".format(eseed))
