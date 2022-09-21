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

err_mq=80
eseed=100

#path="/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_100seeeds_corrhplus/test_seed{1}".format(err_mq,eseed)

#path="/home/td271008/work/cpymadx/Orbit_cpymad/mb_fielderr_roll_mq_offset_{0}_IP5_2it_100seeeds_corrhplus/test_seed{1}".format(err_mq,eseed)

#path="/home/td271008/work/cpymadx/Orbit_cpymad/test_seed{0}".format(eseed)

path="/home/td271008/work/cpymadx/Orbit_cpymad/tune_match_mb_fielderr_roll_mq_offset_{0}_IP5_100seeeds_corrhplus/test_seed{1}".format(err_mq,eseed)


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
fig1.savefig(path+"/orbit_mqoffset_mberrfield_{0}_line_err_seed{1}.pdf".format(err_mq,eseed))

fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(optics_tre['S']/1000., optics_tre['X'], "-r")
ax2[0].set_ylabel("x [m]")
ax2[1].plot(optics_tre['S']/1000., optics_tre['Y'], "-b")
ax2[1].set_xlabel("longitundinal position [km]")
ax2[1].set_ylabel("y [m]")
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
fig2.savefig(path+"/orbit_mqoffset_mberrfield_{0}_line_err_treading_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(eseed)):
    fnameall1=path+"/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall1, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all1=pd.read_csv(fnameall1, skiprows=52, sep='\s+', names=head_opt)
    optics_all1=optics_all1.reset_index(drop=True)
    fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax3[0].plot(optics_all1['S']/1000., optics_all1['X'], "-r")
    ax3[0].set_ylabel("x [m]")
    ax3[1].plot(optics_all1['S']/1000., optics_all1['Y'], "-b")
    ax3[1].set_xlabel("longitundinal position [km]")
    ax3[1].set_ylabel("y [m]")
    fig3, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig3.savefig(path+"/orbit_mqoffset_mberrfield_{0}_line_err_treading_correctall_it1_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(eseed)):
    fnameall2=path+"/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs".format(eseed) 
    head_opt=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all2=pd.read_csv(fnameall2, skiprows=52, sep='\s+', names=head_opt)
    optics_all2=optics_all2.reset_index(drop=True)
    fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax4[0].plot(optics_all2['S']/1000., optics_all2['X'], "-r")
    ax4[0].set_ylabel("x [m]")
    ax4[1].plot(optics_all2['S']/1000., optics_all2['Y'], "-b")
    ax4[1].set_xlabel("longitundinal position [km]")
    ax4[1].set_ylabel("y [m]")
    fig4, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig4.savefig(path+"/orbit_mqoffset_mberrfield_{0}_line_err_treading_correctall_it2_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_tune_match_seed{0}.tfs".format(eseed)):
    fnametune=path+"/FCCee_heb_modett_tune_match_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall2, header=50, sep='\s+', nrows=0).columns[1:]
    optics_tune=pd.read_csv(fnametune, skiprows=52, sep='\s+', names=head_opt)
    optics_tune=optics_tune.reset_index(drop=True)
    fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax5[0].plot(optics_tune['S']/1000., optics_tune['X'], "-r")
    ax5[0].set_ylabel("x [m]")
    ax5[1].plot(optics_tune['S']/1000., optics_tune['Y'], "-b")
    ax5[1].set_xlabel("longitundinal position [km]")
    ax5[1].set_ylabel("y [m]")
    fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
    fig5.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs".format(eseed)):
    fnameall3=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall3, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all3=pd.read_csv(fnameall3, skiprows=52, sep='\s+', names=head_opt)
    optics_all3=optics_all3.reset_index(drop=True)
    fig6, ax6=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax6[0].plot(optics_all3['S']/1000., optics_all3['X'], "-r")
    ax6[0].set_ylabel("x [m]")
    ax6[1].plot(optics_all3['S']/1000., optics_all3['Y'], "-b")
    ax6[1].set_xlabel("longitundinal position [km]")
    ax6[1].set_ylabel("y [m]")
    fig6, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig6.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it1_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it2_seed{0}.tfs".format(eseed)):
    fnameall4=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it2_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall4, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all4=pd.read_csv(fnameall4, skiprows=52, sep='\s+', names=head_opt)
    optics_all4=optics_all4.reset_index(drop=True)
    fig7, ax7=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax7[0].plot(optics_all4['S']/1000., optics_all4['X'], "-r")
    ax7[0].set_ylabel("x [m]")
    ax7[1].plot(optics_all4['S']/1000., optics_all4['Y'], "-b")
    ax7[1].set_xlabel("longitundinal position [km]")
    ax7[1].set_ylabel("y [m]")
    fig7, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig7.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it2_seed{1}.pdf".format(err_mq,eseed))
    
if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it3_seed{0}.tfs".format(eseed)):
    fnameall5=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it3_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall5, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all5=pd.read_csv(fnameall5, skiprows=52, sep='\s+', names=head_opt)
    optics_all5=optics_all5.reset_index(drop=True)
    fig8, ax8=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax8[0].plot(optics_all5['S']/1000., optics_all5['X'], "-r")
    ax8[0].set_ylabel("x [m]")
    ax8[1].plot(optics_all5['S']/1000., optics_all5['Y'], "-b")
    ax8[1].set_xlabel("longitundinal position [km]")
    ax8[1].set_ylabel("y [m]")
    fig8, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig8.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it3_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it4_seed{0}.tfs".format(eseed)):
    fnameall6=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it4_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall6, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all6=pd.read_csv(fnameall6, skiprows=52, sep='\s+', names=head_opt)
    optics_all6=optics_all6.reset_index(drop=True)
    fig9, ax9=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax9[0].plot(optics_all6['S']/1000., optics_all6['X'], "-r")
    ax9[0].set_ylabel("x [m]")
    ax9[1].plot(optics_all6['S']/1000., optics_all6['Y'], "-b")
    ax9[1].set_xlabel("longitundinal position [km]")
    ax9[1].set_ylabel("y [m]")
    fig9, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig9.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it4_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it5_seed{0}.tfs".format(eseed)):
    fnameall7=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it5_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall7, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all7=pd.read_csv(fnameall7, skiprows=52, sep='\s+', names=head_opt)
    optics_all7=optics_all7.reset_index(drop=True)
    fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax10[0].plot(optics_all7['S']/1000., optics_all7['X'], "-r")
    ax10[0].set_ylabel("x [m]")
    ax10[1].plot(optics_all7['S']/1000., optics_all7['Y'], "-b")
    ax10[1].set_xlabel("longitundinal position [km]")
    ax10[1].set_ylabel("y [m]")
    fig10, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    plt.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it5_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it6_seed{0}.tfs".format(eseed)):
    fnameall8=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it6_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall8, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all8=pd.read_csv(fnameall8, skiprows=52, sep='\s+', names=head_opt)
    optics_all8=optics_all8.reset_index(drop=True)
    fig11, ax11=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax11[0].plot(optics_all8['S']/1000., optics_all8['X'], "-r")
    ax11[0].set_ylabel("x [m]")
    ax11[1].plot(optics_all8['S']/1000., optics_all8['Y'], "-b")
    ax11[1].set_xlabel("longitundinal position [km]")
    ax11[1].set_ylabel("y [m]")
    fig11, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig11.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it6_seed{1}.pdf".format(err_mq,eseed))
 
if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it7_seed{0}.tfs".format(eseed)):
    fnameall9=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it7_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall9, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all9=pd.read_csv(fnameall9, skiprows=52, sep='\s+', names=head_opt)
    optics_all9=optics_all9.reset_index(drop=True)
    fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax12[0].plot(optics_all9['S']/1000., optics_all9['X'], "-r")
    ax12[0].set_ylabel("x [m]")
    ax12[1].plot(optics_all9['S']/1000., optics_all9['Y'], "-b")
    ax12[1].set_xlabel("longitundinal position [km]")
    ax12[1].set_ylabel("y [m]")
    fig12, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig12.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it7_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it8_seed{0}.tfs".format(eseed)):
    fnameall10=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it8_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall10, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all10=pd.read_csv(fnameall10, skiprows=52, sep='\s+', names=head_opt)
    optics_all10=optics_all10.reset_index(drop=True)
    fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax13[0].plot(optics_all10['S']/1000., optics_all10['X'], "-r")
    ax13[0].set_ylabel("x [m]")
    ax13[1].plot(optics_all10['S']/1000., optics_all10['Y'], "-b")
    ax13[1].set_xlabel("longitundinal position [km]")
    ax13[1].set_ylabel("y [m]")
    fig13, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig13.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it8_seed{1}.pdf".format(err_mq,eseed))

    
if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it9_seed{0}.tfs".format(eseed)):
    fnameall11=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it9_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall11, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all11=pd.read_csv(fnameall11, skiprows=52, sep='\s+', names=head_opt)
    optics_all11=optics_all11.reset_index(drop=True)
    fig14, ax14=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax14[0].plot(optics_all11['S']/1000., optics_all11['X'], "-r")
    ax14[0].set_ylabel("x [m]")
    ax14[1].plot(optics_all11['S']/1000., optics_all11['Y'], "-b")
    ax14[1].set_xlabel("longitundinal position [km]")
    ax14[1].set_ylabel("y [m]")
    fig14, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig14.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it9_seed{1}.pdf".format(err_mq,eseed))

    
if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it10_seed{0}.tfs".format(eseed)):
    fnameall12=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it10_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall12, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all12=pd.read_csv(fnameall12, skiprows=52, sep='\s+', names=head_opt)
    optics_all12=optics_all12.reset_index(drop=True)
    fig15, ax15=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax15[0].plot(optics_all12['S']/1000., optics_all12['X'], "-r")
    ax15[0].set_ylabel("x [m]")
    ax15[1].plot(optics_all12['S']/1000., optics_all12['Y'], "-b")
    ax15[1].set_xlabel("longitundinal position [km]")
    ax15[1].set_ylabel("y [m]")
    fig15, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig15.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it10_seed{1}.pdf".format(err_mq,eseed))


if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it11_seed{0}.tfs".format(eseed)):
    fnameall13=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it11_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall13, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all13=pd.read_csv(fnameall13, skiprows=52, sep='\s+', names=head_opt)
    optics_all13=optics_all13.reset_index(drop=True)
    fig16, ax16=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax16[0].plot(optics_all13['S']/1000., optics_all13['X'], "-r")
    ax16[0].set_ylabel("x [m]")
    ax16[1].plot(optics_all13['S']/1000., optics_all13['Y'], "-b")
    ax16[1].set_xlabel("longitundinal position [km]")
    ax16[1].set_ylabel("y [m]")
    fig16, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig16.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it11_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it12_seed{0}.tfs".format(eseed)):
    fnameall14=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it12_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall14, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all14=pd.read_csv(fnameall14, skiprows=52, sep='\s+', names=head_opt)
    optics_all14=optics_all14.reset_index(drop=True)
    fig17, ax17=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax17[0].plot(optics_all14['S']/1000., optics_all14['X'], "-r")
    ax17[0].set_ylabel("x [m]")
    ax17[1].plot(optics_all14['S']/1000., optics_all14['Y'], "-b")
    ax17[1].set_xlabel("longitundinal position [km]")
    ax17[1].set_ylabel("y [m]")
    fig17, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig17.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it12_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it13_seed{0}.tfs".format(eseed)):
    fnameall15=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it13_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall15, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all15=pd.read_csv(fnameall15, skiprows=52, sep='\s+', names=head_opt)
    optics_all15=optics_all15.reset_index(drop=True)
    fig18, ax18=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax18[0].plot(optics_all15['S']/1000., optics_all15['X'], "-r")
    ax18[0].set_ylabel("x [m]")
    ax18[1].plot(optics_all15['S']/1000., optics_all15['Y'], "-b")
    ax18[1].set_xlabel("longitundinal position [km]")
    ax18[1].set_ylabel("y [m]")
    fig18, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig18.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it13_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it14_seed{0}.tfs".format(eseed)):
    fnameall16=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it14_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall16, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all16=pd.read_csv(fnameall16, skiprows=52, sep='\s+', names=head_opt)
    optics_all16=optics_all16.reset_index(drop=True)
    fig19, ax19=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax19[0].plot(optics_all16['S']/1000., optics_all16['X'], "-r")
    ax19[0].set_ylabel("x [m]")
    ax19[1].plot(optics_all16['S']/1000., optics_all16['Y'], "-b")
    ax19[1].set_xlabel("longitundinal position [km]")
    ax19[1].set_ylabel("y [m]")
    fig19, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig19.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it14_seed{1}.pdf".format(err_mq,eseed))

if os.path.exists(path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it15_seed{0}.tfs".format(eseed)):
    fnameall17=path+"/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it15_seed{0}.tfs".format(eseed)
    head_opt=pd.read_csv(fnameall17, header=50, sep='\s+', nrows=0).columns[1:]
    optics_all17=pd.read_csv(fnameall17, skiprows=52, sep='\s+', names=head_opt)
    optics_all17=optics_all17.reset_index(drop=True)
    fig20, ax20=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax20[0].plot(optics_all17['S']/1000., optics_all17['X'], "-r")
    ax20[0].set_ylabel("x [m]")
    ax20[1].plot(optics_all17['S']/1000., optics_all17['Y'], "-b")
    ax20[1].set_xlabel("longitundinal position [km]")
    ax20[1].set_ylabel("y [m]")
    fig20, plt.subplots_adjust(left=.20, right=.97, top=.94, bottom=.11)
    fig20.savefig(path+"/orbit_mqoffset_mberrfield_{0}_tune_match_correctall_it15_seed{1}.pdf".format(err_mq,eseed))



