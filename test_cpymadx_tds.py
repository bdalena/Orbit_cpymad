#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
from scipy.stats import norm
from scipy.fftpack import fft
from scipy.fftpack import fftfreq
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os


madx = Madx()

madx.option(echo=True)

os.system("mkdir -p temp")
os.system("mkdir -p results")

madx.set(format="12d")
madx.set(format="20.12g")
madx.set(format="-18s")

madx.option(rbarc=False)

#----- SEQUENCE DEFINITION -----
madx.call(file='optics/FCCee_heb_modett.ele')
madx.call(file='optics/FCCee_heb_modett.seq')
madx.call(file='optics/FCCee_heb_modett.str')

#----- MACRO DEFINITION -----
madx.call(file='FCCee_heb_misc_macros.madx')
madx.call(file='errors/FCCee_errors_macros.madx')

eseed=10
beam_mode=3
injection_mode=1
with_radiate=0
madx.exec('load_beam()')
madx.exec('tws_select()')
tol_tar=1e-12
switch_rf_on=0
vrf400=60  #RF Voltage at 20 GeV
vrf800=0
narc = 8
nsec = 8

os.system("mkdir -p test_seed{0}".format(eseed))

#----- REFERENCE OPTICS -----
madx.command.beam(sequence='fcc_heb', particle='electron')
madx.use(sequence='fcc_heb')
madx.seqedit(sequence='fcc_heb')
madx.flatten()
madx.cycle(start='ip5')
madx.endedit()
madx.use(sequence='fcc_heb')
madx.savebeta(label='BETAIP',place='IP5')
madx.twiss(sequence='fcc_heb')
madx.value('BETAIP->betx')
madx.value('BETAIP->bety')
madx.value('BETAIP->dx')
madx.value('BETAIP->dy')
madx.value('BETAIP->alpx')
madx.value('BETAIP->alpy')
#sys.exit()
#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett.tfs')

#twiss_ref=madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett.tfs')
#T_ref=df.save_twiss(twiss_ref)

madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed))

df.macro_bpm('test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed),"test_seed{0}/redef_bpm.madx".format(eseed),narc,nsec)
madx.call(file='test_seed{0}/redef_bpm.madx'.format(eseed))

df.macro_corr('test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed),"test_seed{0}/redef_corr.madx".format(eseed),narc,nsec)
madx.call(file='test_seed{0}/redef_corr.madx'.format(eseed))

madx.seqedit(sequence='fcc_heb')
madx.flatten()

df.replace_bpm(madx,'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed),narc,nsec)
df.replace_corr(madx,'test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed),narc,nsec)

madx.use(sequence='fcc_heb')

madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_hv_seed{0}.tfs'.format(eseed))

#sys.exit()
'''
MQ_A=df.get_elem(madx.elements,'mq','a')
CORR_A=df.get_elem(madx.elements,'corr','a')
BPM_A=df.get_elem(madx.elements,'bpm','a')
'''

summ=madx.table.summ
qx_tar=summ.q1
qy_tar=summ.q2
dqx_tar=0
dqy_tar=0
          
#----- ERRORS -----
madx.call(file='select_errors.madx')

#----- SAVE ERRORS ----
madx.option(echo=True)
madx.select(flag='error')
madx.esave(file='test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(eseed))

madx.twiss(sequence='fcc_heb', BETA0='BETAIP', file='test_seed{0}/FCCee_heb_modett_err_seed{0}.tfs'.format(eseed)) #line
#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_err.tfs') #ring

#madx.exec('SEXTUOFF')
#madx.twiss(sequence='fcc_heb', BETA0='BETAIP1', file='FCCee_heb_modett_sextoff.tfs') #line
#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_sextoff.tfs') #ring

#----- ANALYTICAL CORRECTORS CALCULATION -----
'''
file1="FCCee_heb_modett.tfs"
file2="FCCee_heb_errors_corr.out"

df.anal_corr_calc(madx,file1,file2)

#----- CORRECT TUNES & CHROMA -----
#twiss_err=madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_err.tfs')
#T_err=df.save_twiss(twiss_err)

print('GLOBAL')
print(madx.globals['k.mcbh.a1.001'])
print(madx.globals['k.mcbh.a1.002'])
print(madx.globals['k.mcbv.a1.001'])
print(madx.globals['k.mcbv.a1.002'])

madx.twiss(sequence='fcc_heb',BETA0='BETAIP1',file='FCCee_heb_modett_cor.tfs') #line
#madx.twiss(sequence='fcc_heb',file='FCCee_heb_modett_cor.tfs') #ring
'''


'''
#tune match
df.tune_match('fcc_heb',1.0E-7)

#chroma match
df.chroma_match('fcc_heb',1.0E-7)
'''

#twiss_cor=madx.twiss(sequence='fcc_heb',file='FCCee_heb_modett_cor.tfs')
#T_cor=df.save_twiss(twiss_cor)
#madx.twiss(sequence='fcc_heb',file='FCCee_heb_modett_cor.tfs')


#madx.twiss(beta0='BETAIP1')

#----- CORRECTION -----

#madx.exec('SEXTUOFF')

tol_err=1e-5


df.svd_first(madx,5,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc5.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,6,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc6.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,7,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc7.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,8,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc8.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,1,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc1.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,2,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc2.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,3,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc3.tfs test_seed{0}/'.format(eseed))
df.svd_first(madx,4,'fcc_heb',0,tol_err,'BETAIP')
os.system('mv FCCee_heb_modett_orbcor_arc4.tfs test_seed{0}/'.format(eseed))


df.svd_ring(madx,'fcc_heb',0,tol_err)
madx.twiss(BETA0='BETAIP', sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_seed{0}.tfs'.format(eseed))


#df.treading_micado(madx,'IP6','IP7','fcc_heb',0,1e-5,'BETAIP1')
#df.treading_micado(madx,'IP7','FCC_HEBIP1_P_','fcc_heb',0,1e-5,'BETAIP1')

#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_orbcor_all.tfs')
#madx.twiss(BETA0='BETAIP1', sequence='fcc_heb', file='FCCee_heb_modett_orbcor_all.tfs')

#twiss3=madx.twiss(beta0='BETAIP1', sequence='fcc_heb', file='FCCee_heb_modett_orbcor.tfs')
#T3=df.save_twiss(twiss3)




#----- PLOT -----

#df.plot_graph(T1,T3)

'''
print('\n')
print(T_ref)
print('\n')

print('\n')
print(T_err)
print('\n')

print('\n')
print(T_cor)
print('\n')
'''

plt.rcParams.update({'font.size':13})
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)


fnameref="test_seed{0}/FCCee_heb_modett_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
optics_ref=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

'''
opt_ref = optics_ref.loc[optics_ref['NAME'].str.contains('BPMN')]
opt_ref = opt_ref.reset_index(drop=True)
print('\n')
print(optics_ref)
print('\n')
'''

fnamerr="test_seed{0}/FCCee_heb_modett_err_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fnamerr, header=50, sep='\s+', nrows=0).columns[1:]
optics_err=pd.read_csv(fnamerr, skiprows=52, sep='\s+', names=head_opt)
optics_err = optics_err.reset_index(drop=True)
print('\n')
print(optics_err)
print('\n')


fname1="test_seed{0}/FCCee_heb_modett_orbcor_arc4.tfs".format(eseed)
head_opt=pd.read_csv(fname1, header=50, sep='\s+', nrows=0).columns[1:]
optics_1=pd.read_csv(fname1, skiprows=52, sep='\s+', names=head_opt)
optics_1 = optics_1.reset_index(drop=True)
print('\n')
print(optics_1)
print('\n')

'''
fname2="FCCee_heb_modett_orbcor_all_sextuoff.tfs"
head_opt=pd.read_csv(fname2, header=50, sep='\s+', nrows=0).columns[1:]
optics_2=pd.read_csv(fname2, skiprows=52, sep='\s+', names=head_opt)
optics_2 = optics_2.reset_index(drop=True)
print('\n')
print(optics_2)
print('\n')
'''

fname3="test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_seed{0}.tfs".format(eseed)
head_opt=pd.read_csv(fname3, header=50, sep='\s+', nrows=0).columns[1:]
optics_3=pd.read_csv(fname3, skiprows=52, sep='\s+', names=head_opt)
optics_3 = optics_3.reset_index(drop=True)
print('\n')
print(optics_3)
print('\n')




#-----------------------------------------------------------------------------
'''
# Plot orbits in m
fig, ax=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax[0].plot(T_err['S']/1000., T_err['X'], "-r")
ax[0].set_ylabel("x [m]")
ax[1].plot(T_err['S']/1000., T_err['Y'], "-b")
ax[1].set_xlabel("longitundinal position [km]")
ax[1].set_ylabel("y [m]")
ax[1].set_xlim(0,93)
fig, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Sans correction T- Orbits in m")

# Plot beta-beat in %
fig1, ax1=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax1[0].plot(T_err['S']/1000., 100.*((T_err['BETX']-T_ref['BETX'])/T_ref['BETX']), "-r")
ax1[0].set_ylabel(r"$\Delta\beta_{x}/\beta_{x,ref}$ [%]")
ax1[1].plot(T_err['S']/1000., 100.*((T_err['BETY']-T_ref['BETY'])/T_ref['BETY']), "-b")
ax1[1].set_xlabel("longitundinal position [km]")
ax1[1].set_ylabel(r"$\Delta\beta_{y}/\beta_{y,ref}$ [%]")
ax1[1].set_xlim(0,93)
fig1, plt.subplots_adjust(left=.13, right=.97, top=.94, bottom=.11)
plt.title("-Sans correction T- Beta-beat in %")

# Plot normalized dispersion mm^1/2
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(T_err['S']/1000.,(T_err['DX']-T_ref['DX'])/np.sqrt(T_ref['BETX']), "-r")
ax2[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
ax2[1].plot(T_err['S']/1000., (T_err['DY']-T_ref['DY'])/np.sqrt(T_ref['BETY']), "-b")
ax2[1].set_xlabel("longitundinal position [km]")
ax2[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
ax2[1].set_xlim(0,93)
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Sans correction T- Normalized dispersion [mm$^{1/2}$]")
'''
#-----------------------------------------------------------------------------
'''
# Plot orbits in m
fig3, ax3=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax3[0].plot(T_cor['S']/1000., T_cor['X'], "-r")
ax3[0].set_ylabel("x [m]")
ax3[1].plot(T_cor['S']/1000., T_cor['Y'], "-b")
ax3[1].set_xlabel("longitundinal position [km]")
ax3[1].set_ylabel("y [m]")
ax3[1].set_xlim(0,93)
fig, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Avec cor analytique T- Orbits in m")

# Plot beta-beat in %
fig4, ax4=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax4[0].plot(T_cor['S']/1000., 100.*((T_cor['BETX']-T_ref['BETX'])/T_ref['BETX']), "-r")
ax4[0].set_ylabel(r"$\Delta\beta_{x}/\beta_{x,ref}$ [%]")
ax4[1].plot(T_cor['S']/1000., 100.*((T_cor['BETY']-T_ref['BETY'])/T_ref['BETY']), "-b")
ax4[1].set_xlabel("longitundinal position [km]")
ax4[1].set_ylabel(r"$\Delta\beta_{y}/\beta_{y,ref}$ [%]")
ax4[1].set_xlim(0,93)
fig4, plt.subplots_adjust(left=.13, right=.97, top=.94, bottom=.11)
plt.title("-Avec cor analytique T- Beta-beat in %")

# Plot normalized dispersion mm^1/2
fig5, ax5=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax5[0].plot(T_cor['S']/1000.,(T_cor['DX']-T_ref['DX'])/np.sqrt(T_ref['BETX']), "-r")
ax5[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
ax5[1].plot(T_cor['S']/1000., (T_cor['DY']-T_ref['DY'])/np.sqrt(T_ref['BETY']), "-b")
ax5[1].set_xlabel("longitundinal position [km]")
ax5[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
ax5[1].set_xlim(0,93)
fig5, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Avec cor analytique T- Normalized dispersion [mm$^{1/2}$]")
'''
#-----------------------------------------------------------------------------

# Plot orbits in m
fig6, ax6=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax6[0].plot(optics_err['S']/1000., optics_err['X'], "-r")
ax6[0].set_ylabel("x [m]")
ax6[1].plot(optics_err['S']/1000., optics_err['Y'], "-b")
ax6[1].set_xlabel("longitundinal position [km]")
ax6[1].set_ylabel("y [m]")
ax6[1].set_xlim(0,93)
fig6, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-with errors (from tfs) (seed={0})- Orbits in m".format(eseed))

'''
# Plot beta-beat in %
fig7, ax7=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax7[0].plot(optics_err['S']/1000., 100.*((optics_err['BETX']-optics_ref['BETX'])/optics_ref['BETX']), "-r")
ax7[0].set_ylabel(r"$\Delta\beta_{x}/\beta_{x,ref}$ [%]")
ax7[1].plot(optics_err['S']/1000., 100.*((optics_err['BETY']-optics_ref['BETY'])/optics_ref['BETY']), "-b")
ax7[1].set_xlabel("longitundinal position [km]")
ax7[1].set_ylabel(r"$\Delta\beta_{y}/\beta_{y,ref}$ [%]")
ax7[1].set_xlim(0,93)
fig7, plt.subplots_adjust(left=.13, right=.97, top=.94, bottom=.11)
plt.title("-Sans correction DF- Beta-beat in %")

# Plot normalized dispersion mm^1/2
fig8, ax8=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax8[0].plot(optics_err['S']/1000.,(optics_err['DX']-optics_ref['DX'])/np.sqrt(optics_ref['BETX']), "-r")
ax8[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
ax8[1].plot(optics_err['S']/1000., (optics_err['DY']-optics_ref['DY'])/np.sqrt(optics_ref['BETY']), "-b")
ax8[1].set_xlabel("longitundinal position [km]")
ax8[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
ax8[1].set_xlim(0,93)
fig8, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Sans correction DF- Normalized dispersion [mm$^{1/2}$]")
'''
#-----------------------------------------------------------------------------
'''
# Plot orbits in m
fig9, ax9=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax9[0].plot(optics_cor['S']/1000., optics_cor['X'], "-r")
ax9[0].set_ylabel("x [m]")
ax9[1].plot(optics_cor['S']/1000., optics_cor['Y'], "-b")
ax9[1].set_xlabel("longitundinal position [km]")
ax9[1].set_ylabel("y [m]")
ax9[1].set_xlim(0,93)
fig9, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-with analytical corr (from tfs)- Orbits in m")

# Plot beta-beat in %
fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax10[0].plot(optics_cor['S']/1000., 100.*((optics_cor['BETX']-optics_ref['BETX'])/optics_ref['BETX']), "-r")
ax10[0].set_ylabel(r"$\Delta\beta_{x}/\beta_{x,ref}$ [%]")
ax10[1].plot(optics_cor['S']/1000., 100.*((optics_cor['BETY']-optics_ref['BETY'])/optics_ref['BETY']), "-b")
ax10[1].set_xlabel("longitundinal position [km]")
ax10[1].set_ylabel(r"$\Delta\beta_{y}/\beta_{y,ref}$ [%]")
ax10[1].set_xlim(0,93)
fig10, plt.subplots_adjust(left=.13, right=.97, top=.94, bottom=.11)
plt.title("-Avec cor analytique DF- Beta-beat in %")

# Plot normalized dispersion mm^1/2
fig11, ax11=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax11[0].plot(optics_cor['S']/1000.,(optics_cor['DX']-optics_ref['DX'])/np.sqrt(optics_ref['BETX']), "-r")
ax11[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
ax11[1].plot(optics_cor['S']/1000., (optics_cor['DY']-optics_ref['DY'])/np.sqrt(optics_ref['BETY']), "-b")
ax11[1].set_xlabel("longitundinal position [km]")
ax11[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
ax11[1].set_xlim(0,93)
fig11, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-Avec cor analytique DF- Normalized dispersion [mm$^{1/2}$]")
'''

#-----------------------------------------------------------------------------
'''
# Plot orbits in m
fig10, ax10=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax10[0].plot(optics_sext['S']/1000., optics_sext['X'], "-r")
ax10[0].set_ylabel("x [m]")
ax10[1].plot(optics_sext['S']/1000., optics_sext['Y'], "-b")
ax10[1].set_xlabel("longitundinal position [km]")
ax10[1].set_ylabel("y [m]")
ax10[1].set_xlim(0,93)
fig10, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-with sextoff (from tfs)- Orbits in m")
'''


# Plot orbits in m
fig12, ax12=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax12[0].plot(optics_1['S']/1000., optics_1['X'], "-r")
ax12[0].set_ylabel("x [m]")
ax12[1].plot(optics_1['S']/1000., optics_1['Y'], "-b")
ax12[1].set_xlabel("longitundinal position [km]")
ax12[1].set_ylabel("y [m]")
ax12[1].set_xlim(0,93)
fig12, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-End of treading (seed={0})- Orbits in m".format(eseed))

'''
fig13, ax13=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax13[0].plot(optics_2['S']/1000., optics_2['X'], "-r")
ax13[0].set_ylabel("x [m]")
ax13[1].plot(optics_2['S']/1000., optics_2['Y'], "-b")
ax13[1].set_xlabel("longitundinal position [km]")
ax13[1].set_ylabel("y [m]")
ax13[1].set_xlim(0,93)
fig13, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-with CORR madx and sextuoff- Orbits in m")
'''

fig14, ax14=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax14[0].plot(optics_3['S']/1000., optics_3['X'], "-r")
ax14[0].set_ylabel("x [m]")
ax14[1].plot(optics_3['S']/1000., optics_3['Y'], "-b")
ax14[1].set_xlabel("longitundinal position [km]")
ax14[1].set_ylabel("y [m]")
ax14[1].set_xlim(0,93)
fig14, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("-with CORR madx and sextuon (seed={0})- Orbits in m".format(eseed))




plt.show()
