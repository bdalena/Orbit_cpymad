#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import def_functions2 as df
import numpy as np
import pandas as pd
import sys
import os

madx = Madx()

madx.option(echo=True)

madx.set(format="12d")
madx.set(format="20.12g")
madx.set(format="-18s")

madx.option(rbarc=False)

#----- SEQUENCE DEFINITION -----
madx.call(file='/local/home/bdalena/FCC-ee/FCC-ee-lattice/check_optics/fcc-ee-heb/optics/4IPs/optics_PA31/merge_orbit_cpymad/Orbit_cpymad/optics/FCCee_heb_modett.ele')
madx.call(file='/local/home/bdalena/FCC-ee/FCC-ee-lattice/check_optics/fcc-ee-heb/optics/4IPs/optics_PA31/merge_orbit_cpymad/Orbit_cpymad/optics/FCCee_heb_modett.seq')
madx.call(file='/local/home/bdalena/FCC-ee/FCC-ee-lattice/check_optics/fcc-ee-heb/optics/4IPs/optics_PA31/merge_orbit_cpymad/Orbit_cpymad/optics/FCCee_heb_modett.str')

#----- MACRO DEFINITION -----
madx.call(file='/local/home/bdalena/FCC-ee/FCC-ee-lattice/check_optics/fcc-ee-heb/optics/4IPs/optics_PA31/merge_orbit_cpymad/Orbit_cpymad/FCCee_heb_misc_macros.madx')
madx.call(file='/local/home/bdalena/FCC-ee/FCC-ee-lattice/check_optics/fcc-ee-heb/optics/4IPs/optics_PA31/merge_orbit_cpymad/Orbit_cpymad/errors/FCCee_errors_macros_girder.madx')

eseed=47
beam_mode=3
injection_mode=0
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

summ=madx.table.summ
qx_tar=summ.q1
qy_tar=summ.q2
dqx_tar=0
dqy_tar=0
          
#----- ERRORS -----
#----- Seed for errors assignment ----
madx.input("EOPTION, seed=47;")
madx.call(file='select_errors_girder.madx')

#----- SAVE ERRORS ----
madx.option(echo=True)
madx.select(flag='error')
madx.esave(file='test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(eseed))

madx.twiss(sequence='fcc_heb', BETA0='BETAIP', file='test_seed{0}/FCCee_heb_modett_err_seed{0}.tfs'.format(eseed)) #line
#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_err.tfs') #ring



#----- ANALYTICAL CORRECTORS CALCULATION -----
'''
file1='test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed)
file2='test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(eseed)

df.anal_corr_calc(file1,file2)
'''

#----- CORRECTION -----

madx.exec('SEXTUOFF')

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

#SVD line: iteration n°1
df.svd_ring(madx,'fcc_heb',0,tol_err)
#madx.twiss(BETA0='BETAIP', sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it1_seed{0}.tfs'.format(eseed))
madx.twiss(BETA0='BETAIP', sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it1_seed{0}.tfs'.format(eseed))

#SVD line: iteration n°2
df.svd_ring(madx,'fcc_heb',0,tol_err)
#madx.twiss(BETA0='BETAIP', sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{0}.tfs'.format(eseed))
madx.twiss(BETA0='BETAIP', sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_line_sextuoff_it2_seed{0}.tfs'.format(eseed))

os.system('mv cx_fccee_heb_mic_all.tab cx_fccee_heb_mic_all_line_sextuoff_it2.tab')
os.system('mv cy_fccee_heb_mic_all.tab cy_fccee_heb_mic_all_line_sextuoff_it2.tab')

#SVD ring: iteration n°1
df.svd_ring(madx,'fcc_heb',0,tol_err)
madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs'.format(eseed))


fname='test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it1_seed{0}.tfs'.format(eseed)
head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
optics=pd.read_csv(fname, skiprows=52, sep='\s+', names=head_opt)
optics=optics.reset_index(drop=True)
rms_orbit_x=np.std(optics['X'])
rms_orbit_y=np.std(optics['Y'])

file1='test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed)
file2='test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(eseed)
mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)

max_it=15
nb_it=1

while rms_orbit_x>orbit_x or rms_orbit_y>orbit_y and nb_it<max_it:
    
    nb_it+=1
    df.svd_ring(madx,'fcc_heb',0,tol_err)
    madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{1}_seed{0}.tfs'.format(eseed,nb_it))
    os.system('mv cx_fccee_heb_mic_all.tab cx_fccee_heb_mic_all_ring_sextuoff_it{0}.tab'.format(nb_it))
    os.system('mv cy_fccee_heb_mic_all.tab cy_fccee_heb_mic_all_ring_sextuoff_it{0}.tab'.format(nb_it))

    fname='test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuoff_it{1}_seed{0}.tfs'.format(eseed,nb_it)
    head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
    optics=pd.read_csv(fname, skiprows=52, sep='\s+', names=head_opt)
    optics=optics.reset_index(drop=True)
    rms_orbit_x=np.std(optics['X'])
    rms_orbit_y=np.std(optics['Y'])


    
madx.exec('SEXTUON')

#SVD ring: iteration n°1
df.svd_ring(madx,'fcc_heb',0,tol_err)
madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{0}.tfs'.format(eseed))


sys.exit()

"""
We stop here since the matching procedure is not working for the moment.
"""

df.match_quad(madx,'test_seed{0}/FCCee_heb-quads.str')

madx.call(file='test_seed{0}/FCCee_heb-quads.str')

madx.command.match(sequence='fcc_heb')
madx.vary(name='k1f_tune',step=1.0E-7)
madx.vary(name='k1d_tune',step=1.0E-7)
madx.command.global_(sequence='fcc_heb',Q1=415.225,Q2=416.29)
madx.jacobian(calls=10,tolerance=tol_tar)
madx.endmatch()

madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_tune_match_seed{0}.tfs'.format(eseed))

    
sys.exit()

#----- MATCHING -----

df.match_quad(madx,'test_seed{0}/FCCee_heb-quads.str')

madx.call(file='test_seed{0}/FCCee_heb-quads.str')

madx.command.match(sequence='fcc_heb')
madx.vary(name='k1f_tune',step=1.0E-7)
madx.vary(name='k1d_tune',step=1.0E-7)
madx.command.global_(sequence='fcc_heb',Q1=415.225,Q2=416.29)
madx.jacobian(calls=10,tolerance=tol_tar)
madx.endmatch()

madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_tune_match_seed{0}.tfs'.format(eseed))

#iteration n°1 = initialisation
df.svd_ring(madx,'fcc_heb',0,tol_err)
madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs'.format(eseed))
os.system('mv cx_fccee_heb_mic_all.tab cx_fccee_heb_mic_all_ring_it1.tab')
os.system('mv cy_fccee_heb_mic_all.tab cy_fccee_heb_mic_all_ring_it1.tab')

fname='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it1_seed{0}.tfs'.format(eseed)
head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
optics=pd.read_csv(fname, skiprows=52, sep='\s+', names=head_opt)
optics=optics.reset_index(drop=True)
rms_orbit_x=np.std(optics['X'])
rms_orbit_y=np.std(optics['Y'])

file1='test_seed{0}/FCCee_heb_modett_seed{0}.tfs'.format(eseed)
file2='test_seed{0}/FCCee_heb_errors_corr_seed{0}.out'.format(eseed)
mean_corr_x,mean_corr_y,orbit_x,orbit_y=df.anal_corr_calc(file1,file2)

max_it=15
nb_it=1

while rms_orbit_x>orbit_x or rms_orbit_y>orbit_y and nb_it<max_it:
    
    nb_it+=1
    df.svd_ring(madx,'fcc_heb',0,tol_err)
    madx.twiss(sequence='fcc_heb', file='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs'.format(eseed,nb_it))
    os.system('mv cx_fccee_heb_mic_all.tab cx_fccee_heb_mic_all_ring_it{0}.tab'.format(nb_it))
    os.system('mv cy_fccee_heb_mic_all.tab cy_fccee_heb_mic_all_ring_it{0}.tab'.format(nb_it))

    fname='test_seed{0}/FCCee_heb_modett_orbcor_all_sextuon_tune_match_it{1}_seed{0}.tfs'.format(eseed,nb_it)
    head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
    optics=pd.read_csv(fname, skiprows=52, sep='\s+', names=head_opt)
    optics=optics.reset_index(drop=True)
    rms_orbit_x=np.std(optics['X'])
    rms_orbit_y=np.std(optics['Y'])
    
