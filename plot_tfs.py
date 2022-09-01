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

err_mq=20
eseed=100

path="/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_100seeeds_corrhplus/test_seed{1}".format(err_mq,eseed)

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

fnameall1=path+"/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt)
optics_all1=optics_all1.reset_index(drop=True)

fnameall2=path+"/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
optics_all2=pd.read_csv(fnameall2, skiprows=52, sep='\s+', names=head_opt)
optics_all2=optics_all2.reset_index(drop=True)

#-----PLOT-----

# Plot orbits in m
fig1, ax1=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax1[0].plot(optics_err['S']/1000., optics_err['X'], "-r")
ax1[0].set_ylabel("x [m]")
ax1[1].plot(optics_err['S']/1000., optics_err['Y'], "-b")
ax1[1].set_xlabel("longitundinal position [km]")
ax1[1].set_ylabel("y [m]")
fig1, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.savefig(path+"/orbit_mqoffset_{0}_line_err_seed{1}.pdf".format(err_mq,eseed))

fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(optics_tre['S']/1000., optics_tre['X'], "-r")
ax2[0].set_ylabel("x [m]")
ax2[1].plot(optics_tre['S']/1000., optics_tre['Y'], "-b")
ax2[1].set_xlabel("longitundinal position [km]")
ax2[1].set_ylabel("y [m]")
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.savefig(path+"/orbit_mqoffset_{0}_line_err_treading_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(fnameall1):
    fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax3[0].plot(optics_all1['S']/1000., optics_all1['X'], "-r")
    ax3[0].set_ylabel("x [m]")
    ax3[1].plot(optics_all1['S']/1000., optics_all1['Y'], "-b")
    ax3[1].set_xlabel("longitundinal position [km]")
    ax3[1].set_ylabel("y [m]")
    fig3, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    plt.savefig(path+"/orbit_mqoffset_{0}_line_err_treading_correctall_it1_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(fnameall2):
    fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax4[0].plot(optics_all2['S']/1000., optics_all2['X'], "-r")
    ax4[0].set_ylabel("x [m]")
    ax4[1].plot(optics_all2['S']/1000., optics_all2['Y'], "-b")
    ax4[1].set_xlabel("longitundinal position [km]")
    ax4[1].set_ylabel("y [m]")
    fig4, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    plt.savefig(path+"/orbit_mqoffset_{0}_line_err_treading_correctall_it2_seed{1}.pdf".format(err_mq,eseed))




