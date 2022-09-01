#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.fftpack import fft   
import numpy as np
import pandas as pd
import sys
import os


#madx = Madx()

#----------------------------------------------------------------------------------------------------------------

def std(data,elem1,elem2):
    moy_x = np.mean(data[elem1])
    std_x = np.std(data[elem1])
    moy_y = np.mean(data[elem2])
    std_y = np.std(data[elem2])
    return(moy_x,moy_y,std_x,std_y)

#- - - - - -

def store_corr(data,ind_k1,strength_x,strength_y):
    corr_x=[]
    corr_y=[]
    for p in range(len(data)):
        if data.iloc[p,ind_k1]>0:
            corr_x.append(str(strength_x))
            corr_y.append(str(0.0))
        else:
            corr_y.append(str(strength_y))
            corr_x.append(str(0.0))
    return(corr_x,corr_y)

#- - - - - -

def global_corr(madx,arc,corr_x,corr_y):
    for i in range(1,370,1):
        if i<10:
            madx.globals['k.mcbh.a{0}.00{1}'.format(arc,i)]=corr_x[i-1]
            madx.globals['k.mcbv.a{0}.00{1}'.format(arc,i)]=corr_y[i-1]
        elif i>=10 and i<100:
            madx.globals['k.mcbh.a{0}.0{1}'.format(arc,i)]=corr_x[i-1]
            madx.globals['k.mcbv.a{0}.0{1}'.format(arc,i)]=corr_y[i-1]
        elif i>=100 and i<370:
            madx.globals['k.mcbh.a{0}.{1}'.format(arc,i)]=corr_x[i-1]
            madx.globals['k.mcbv.a{0}.{1}'.format(arc,i)]=corr_y[i-1]

#- - - - - -

def anal_corr_calc(madx,file1,file2):
    '''
    Gives a first (analytical) value for the correctors strength.
    Generates a str file with with the calculated correctors strength.
    madx = the object madx
    file1 = .tfs of your reference 
    file2 = .tfs of your result with errors
    '''
    
    plt.rcParams.update({'font.size':13})
    plt.rc('xtick',labelsize=11)
    plt.rc('ytick',labelsize=11)

    fnameref1=file1
    head_opt1=pd.read_csv(fnameref1, header=50, sep='\s+', nrows=0).columns[1:]
    ref_tfs=pd.read_csv(fnameref1, skiprows=52, sep='\s+', names=head_opt1)
    ref_tfs_headers=pd.read_csv(fnameref1, sep='\s+', names=head_opt1)

    fnameref2=file2
    head_opt2=pd.read_csv(fnameref2, header=6, sep='\s+', nrows=0).columns[1:]
    ref_err=pd.read_csv(fnameref2, skiprows=8, sep='\s+', names=head_opt2)

    s_cell = "MS.A1.021"
    e_cell = "MS.A1.023"
    ind1 = ref_tfs.index[ref_tfs["NAME"]==s_cell].tolist()[0]
    ind2 = ref_tfs.index[ref_tfs["NAME"]==e_cell].tolist()[0]

    #----- VARIABLES -----

    ind_q=ref_tfs.index[ref_tfs["NAME"]=="MQ.A1.021"].tolist()[0] #index of the first quad selected in the arc
    ind_d=ref_tfs.index[ref_tfs["NAME"]=="MB.A1.021A"].tolist()[0] #index of the first dpole selected in the arc
    ind_S_arc=ref_tfs.index[ref_tfs["NAME"]=="MQ.A1.001"].tolist()[0] #index for the begging of an arc
    ind_E_arc=ref_tfs.index[ref_tfs["NAME"]=="MQ.A2.001"].tolist()[0] #index for the ending of an arc

    #Getting back the std
    data_MQ_A1_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A1')]
    data_MQ_A2_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A2')]
    data_MQ_A3_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A3')]
    data_MQ_A4_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A4')]
    data_MQ_A5_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A5')]
    data_MQ_A6_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A6')]
    data_MQ_A7_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A7')]
    data_MQ_A8_err=ref_err.loc[ref_err["NAME"].str.contains('MQ.A8')]

    data_MB_A1_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A1')]
    data_MB_A2_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A2')]
    data_MB_A3_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A3')]
    data_MB_A4_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A4')]
    data_MB_A5_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A5')]
    data_MB_A6_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A6')]
    data_MB_A7_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A7')]
    data_MB_A8_err=ref_err.loc[ref_err["NAME"].str.contains('MB.A8')]
    
    data_corr=ref_tfs.loc[ref_tfs["NAME"].str.contains('CORR.A')]
    data_bpm=ref_tfs.loc[ref_tfs["NAME"].str.contains('BPM')]
    data_MQ=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ')]
    data_MQA=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A')]
    data_MQS=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.S')]
    data_MB=ref_tfs.loc[ref_tfs["NAME"].str.contains('MB')]
    data_MS=ref_tfs.loc[ref_tfs["NAME"].str.contains('MS')]

    ma1x,ma1y,std1x,std1y=std(data_MQ_A1_err,"DX","DY")
    ma2x,ma2y,std2x,std2y=std(data_MQ_A2_err,"DX","DY")
    ma3x,ma3y,std3x,std3y=std(data_MQ_A3_err,"DX","DY")
    ma4x,ma4y,std4x,std4y=std(data_MQ_A4_err,"DX","DY")
    ma5x,ma5y,std5x,std5y=std(data_MQ_A5_err,"DX","DY")
    ma6x,ma6y,std6x,std6y=std(data_MQ_A6_err,"DX","DY")
    ma7x,ma7y,std7x,std7y=std(data_MQ_A7_err,"DX","DY")
    ma8x,ma8y,std8x,std8y=std(data_MQ_A8_err,"DX","DY")

    ma1x,ma1y,std1x_mb,std1y_mb=std(data_MB_A1_err,"K0L","DPSI")
    ma2x,ma2y,std2x_mb,std2y_mb=std(data_MB_A2_err,"K0L","DPSI")
    ma3x,ma3y,std3x_mb,std3y_mb=std(data_MB_A3_err,"K0L","DPSI")
    ma4x,ma4y,std4x_mb,std4y_mb=std(data_MB_A4_err,"K0L","DPSI")
    ma5x,ma5y,std5x_mb,std5y_mb=std(data_MB_A5_err,"K0L","DPSI")
    ma6x,ma6y,std6x_mb,std6y_mb=std(data_MB_A6_err,"K0L","DPSI")
    ma7x,ma7y,std7x_mb,std7y_mb=std(data_MB_A7_err,"K0L","DPSI")
    ma8x,ma8y,std8x_mb,std8y_mb=std(data_MB_A8_err,"K0L","DPSI")

    #Columns index
    i_S=ref_tfs.columns.get_loc("S")
    i_L=ref_tfs.columns.get_loc("L")
    i_K1L=ref_tfs.columns.get_loc("K1L")
    i_ANGLE=ref_tfs.columns.get_loc("ANGLE")
    i_DPSI=ref_err.columns.get_loc("DPSI")

    L_tot=ref_tfs.iat[len(ref_tfs)-1,i_S]-ref_tfs.iat[0,i_S] # (m)
    L_cell=ref_tfs.iat[ind2,i_S]-ref_tfs.iat[ind1,i_S] # (m)
    L_d=ref_tfs.iat[ind_d,i_L] #dipole length (m)
    L_q=ref_tfs.iat[ind_q,i_L] #quad length (m)
    K1L=ref_tfs.iat[ind_q,i_K1L]
    MB_angle=ref_tfs.iat[ind_d,i_ANGLE] #dipole bending angle (rad)

    #from headers
    p=float(ref_tfs_headers.iat[7,3])
    Q_x=float(ref_tfs_headers.iat[22,3])
    Q_y=float(ref_tfs_headers.iat[23,3])

    Brho=3.3356*p
    n=4 #number of dipoles in the cell

    MB_field_err_a1=std1x_mb/L_d
    MB_field_err_a2=std2x_mb/L_d
    MB_field_err_a3=std3x_mb/L_d
    MB_field_err_a4=std4x_mb/L_d
    MB_field_err_a5=std5x_mb/L_d
    MB_field_err_a6=std6x_mb/L_d
    MB_field_err_a7=std7x_mb/L_d
    MB_field_err_a8=std8x_mb/L_d

    mu_cell=np.pi/2
    beta_bar=L_cell/np.sin(mu_cell)
    beta_max=beta_bar*(1+np.sin(L_cell/2))
    rho=L_tot/(2*np.pi) #radius of the circle
    N_d=(2*np.pi)/MB_angle #tot nb of dipole
    N_q=N_d/2 #tot nb of quad (because we work in the arcs)

    #----- CORRECTOR STRENGTH -----
    
    f=1
    fact1x=f#*np.sign(ma1x)
    fact2x=f#*np.sign(ma2x)
    fact3x=f#*np.sign(ma3x)
    fact4x=f#*np.sign(ma4x)
    fact5x=f#*np.sign(ma5x)
    fact6x=f#*np.sign(ma6x)
    fact7x=f#*np.sign(ma7x)
    fact8x=f#*np.sign(ma8x)

    #Formula 1: from report 43
    corr_strength_x_A1=fact1x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a1)**2+2*std1x**2*K1L**2))
    corr_strength_x_A2=fact2x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a2)**2+2*std2x**2*K1L**2))
    corr_strength_x_A3=fact3x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a3)**2+2*std3x**2*K1L**2))
    corr_strength_x_A4=fact4x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a4)**2+2*std4x**2*K1L**2))
    corr_strength_x_A5=fact5x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a5)**2+2*std5x**2*K1L**2))
    corr_strength_x_A6=fact6x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a6)**2+2*std6x**2*K1L**2))
    corr_strength_x_A7=fact7x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a7)**2+2*std7x**2*K1L**2))
    corr_strength_x_A8=fact8x*np.sqrt((beta_bar/beta_max)*(n*(MB_angle*MB_field_err_a8)**2+2*std8x**2*K1L**2))
    
    fact1y=f#*np.sign(ma1y)
    fact2y=f#*np.sign(ma2y)
    fact3y=f#*np.sign(ma3y)
    fact4y=f#*np.sign(ma4y)
    fact5y=f#*np.sign(ma5y)
    fact6y=f#*np.sign(ma6y)
    fact7y=f#*np.sign(ma7y)
    fact8y=f#*np.sign(ma8y)

    #Formula 1: from report 43
    corr_strength_y_A1=fact1y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std1y_mb/rho)**2+2*std1y**2*K1L**2))
    corr_strength_y_A2=fact2y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std2y_mb/rho)**2+2*std2y**2*K1L**2))
    corr_strength_y_A3=fact3y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std3y_mb/rho)**2+2*std3y**2*K1L**2))
    corr_strength_y_A4=fact4y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std4y_mb/rho)**2+2*std4y**2*K1L**2))
    corr_strength_y_A5=fact5y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std5y_mb/rho)**2+2*std5y**2*K1L**2))
    corr_strength_y_A6=fact6y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std6y_mb/rho)**2+2*std6y**2*K1L**2))
    corr_strength_y_A7=fact7y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std7y_mb/rho)**2+2*std7y**2*K1L**2))
    corr_strength_y_A8=fact8y*np.sqrt((beta_bar/beta_max)*(n*(L_d*std8y_mb/rho)**2+2*std8y**2*K1L**2))
    
    #----- STORE CORR VALUE -----

    data_MQ_A1=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A1')]
    data_MQ_A2=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A2')]
    data_MQ_A3=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A3')]
    data_MQ_A4=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A4')]
    data_MQ_A5=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A5')]
    data_MQ_A6=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A6')]
    data_MQ_A7=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A7')]
    data_MQ_A8=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A8')]
  
    corr_x_A1,corr_y_A1=store_corr(data_MQ_A1,i_K1L,corr_strength_x_A1,corr_strength_y_A1)
    corr_x_A2,corr_y_A2=store_corr(data_MQ_A2,i_K1L,corr_strength_x_A2,corr_strength_y_A2)
    corr_x_A3,corr_y_A3=store_corr(data_MQ_A3,i_K1L,corr_strength_x_A3,corr_strength_y_A3)
    corr_x_A4,corr_y_A4=store_corr(data_MQ_A4,i_K1L,corr_strength_x_A4,corr_strength_y_A4)
    corr_x_A5,corr_y_A5=store_corr(data_MQ_A5,i_K1L,corr_strength_x_A5,corr_strength_y_A5)
    corr_x_A6,corr_y_A6=store_corr(data_MQ_A6,i_K1L,corr_strength_x_A6,corr_strength_y_A6)
    corr_x_A7,corr_y_A7=store_corr(data_MQ_A7,i_K1L,corr_strength_x_A7,corr_strength_y_A7)
    corr_x_A8,corr_y_A8=store_corr(data_MQ_A8,i_K1L,corr_strength_x_A8,corr_strength_y_A8)
  
    global_corr(madx,1,corr_x_A1,corr_y_A1)
    global_corr(madx,2,corr_x_A2,corr_y_A2)
    global_corr(madx,3,corr_x_A3,corr_y_A3)
    global_corr(madx,4,corr_x_A4,corr_y_A4)
    global_corr(madx,5,corr_x_A5,corr_y_A5)
    global_corr(madx,6,corr_x_A6,corr_y_A6)
    global_corr(madx,7,corr_x_A7,corr_y_A7)
    global_corr(madx,8,corr_x_A8,corr_y_A8)

    return

#----------------------------------------------------------------------------------------------------------------  

def get_elem(madx_elem,elem_name,elem_loc):
    '''
    Saving the characteristics of an element type.
    madx_elem = madx.elements (all elements with their characteristcs)
    elem_name = name of the element : str ('mq','corr','bpm'...)
    elem_loc = localisation of the element : str (arc:'a' or section:'s')
    return a string with the name of the elements and their characteristics.
    '''
    raw_list=[]
    for j in range(1,9):
        i=1
        while i<10:
            raw_list.append(madx_elem.get(elem_name+"."+elem_loc+"{0}.00{1}".format(j,i)))
            i+=1

        i=10
        while i<100:
            raw_list.append(madx_elem.get(elem_name+"."+elem_loc+"{0}.0{1}".format(j,i)))
            i+=1

        i=100
        while i<370:
            raw_list.append(madx_elem.get(elem_name+"."+elem_loc+"{0}.{1}".format(j,i)))
            i+=1
    final_str='\n'.join(map(str,raw_list))
    return(final_str)

#----------------------------------------------------------------------------------------------------------------  

def tune_match(seq,pace):
    '''
    For matching the tune.
    seq = name of the sequence : str
    pace = step for varying : float
    '''
    madx.command.match(sequence=seq)
    madx.vary(name='k1f_tune',step=pace)
    madx.vary(name='k1d_tune',step=pace)
    #madx.command.global_(sequence='fcc_heb',q1=qx_tar,q2=qy_tar) #!!!
    madx.jacobian(calls=10,tolerance=tol_tar)
    madx.endmatch()

#----------------------------------------------------------------------------------------------------------------  

def chroma_match(seq,pace):
    '''
    For matching the chroma.
    seq = name of the sequence : str
    pace = step for varying : float
    '''
    madx.command.match(sequence=seq)
    madx.vary(name='k2sf',step=pace)
    madx.vary(name='k2sd',step=pace)
    #madx.command.global_(sequence='fcc_heb',dq1=dqx_tar,dq2=dqy_tar) #!!!
    madx.jacobian(calls=10,tolerance=tol_tar)
    madx.endmatch()

#----------------------------------------------------------------------------------------------------------------

def save_twiss(twiss):
    twiss_s=twiss.s
    twiss_x=twiss.x
    twiss_y=twiss.y
    twiss_betx=twiss.betx
    twiss_bety=twiss.bety
    twiss_dx=twiss.dx
    twiss_dy=twiss.dy

    T = pd.DataFrame(list(zip(twiss_s,twiss_x,twiss_y,twiss.betx,twiss.bety,twiss.dx,twiss.dx)), columns = ['S','X','Y','BETX','BETY','DX','DY'])

    return(T)

#----------------------------------------------------------------------------------------------------------------
def rename_bpm_sec(data,ind_k1,file_name,nsec):
    for p in range(len(data)):
        jj=p%len(data)+1
       # if jj==1:
       #     ii+=1

        if data.iloc[p,ind_k1]>0:
            if jj<10:
                file_name.write("bpmn.s{0}.00{1} : hbpm;\n".format(nsec,jj))
            elif jj>=10 and jj<100:
                file_name.write("bpmn.s{0}.0{1} : hbpm;\n".format(nsec,jj))

        elif data.iloc[p,ind_k1]<0:
            if jj<10:
                file_name.write("bpmn.s{0}.00{1} : vbpm;\n".format(nsec,jj))
            elif jj>=10 and jj<100:
                file_name.write("bpmn.s{0}.0{1} : vbpm;\n".format(nsec,jj))
                   
def rename_bpm_arc(data,ind_k1,file_name,narc):
    ii=0
    for p in range(len(data)):
        jj=int(p%(len(data)/narc)+1)
        if jj==1:
            ii+=1
            
        if data.iloc[p,ind_k1]>0:
            if jj<10:
                file_name.write("bpmn.a{0}.00{1} : hbpm;\n".format(ii,jj))
            elif jj>=10 and jj<100:
                file_name.write("bpmn.a{0}.0{1} : hbpm;\n".format(ii,jj))
            elif jj>=100 and jj<370:
                file_name.write("bpmn.a{0}.{1} : hbpm;\n".format(ii,jj))
                
        elif data.iloc[p,ind_k1]<0:
            if jj<10:
                file_name.write("bpmn.a{0}.00{1} : vbpm;\n".format(ii,jj))
            elif jj>=10 and jj<100:
                file_name.write("bpmn.a{0}.0{1} : vbpm;\n".format(ii,jj))
            elif jj>=100 and jj<370:
                file_name.write("bpmn.a{0}.{1} : vbpm;\n".format(ii,jj))

#- - - - - -
                
def replace_bpm(madx,file_ref,narc,nsec):

    fnameref=file_ref
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    ref_tfs=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

    data_MQ_A=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A')]
    
    ii=0
    for p in range(len(data_MQ_A)):
        jj=int(p%(len(data_MQ_A)/narc)+1)
        if jj==1:
            ii+=1

        if jj<10:
            madx.replace(element='bpm.a{0}.00{1}'.format(ii,jj), by='bpmn.a{0}.00{1}'.format(ii,jj))
        elif jj>=10 and jj<100:
            madx.replace(element='bpm.a{0}.0{1}'.format(ii,jj), by='bpmn.a{0}.0{1}'.format(ii,jj))
        elif jj>=100 and jj<370:
            madx.replace(element='bpm.a{0}.{1}'.format(ii,jj), by='bpmn.a{0}.{1}'.format(ii,jj))

    for ii in range(0,nsec,1): 
        data_MQ_S=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.S{0}'.format(ii+1))]

        #ii=0
        for p in range(len(data_MQ_S)):
            jj=p%len(data_MQ_S)+1
            if jj==1:
                ii+=1

            if jj<10:
                madx.replace(element='bpm.s{0}.00{1}'.format(ii,jj), by='bpmn.s{0}.00{1}'.format(ii,jj))
            elif jj>=10 and jj<100:
                madx.replace(element='bpm.s{0}.0{1}'.format(ii,jj), by='bpmn.s{0}.0{1}'.format(ii,jj))
            elif jj>=100 and jj<370:
                madx.replace(element='bpm.s{0}.{1}'.format(ii,jj), by='bpmn.s{0}.{1}'.format(ii,jj))
   
    
#- - - - - -

def macro_bpm(file_ref,file_out,narc,nsec):

    fnameref=file_ref
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    ref_tfs=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

    if os.path.exists(file_out):
        os.remove(file_out)
    else:
        print('Not possible: the file does not exist')

    fname=open(file_out, "w+")
    
    #fname.write("HVBPM : macro ={\n")

    fname.write("hbpm : hmonitor, l:=l_bpm;\n")
    fname.write("vbpm : vmonitor, l:=l_bpm;\n")

    data_MQ_A=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A')]
    ia_K1L=ref_tfs.columns.get_loc("K1L")
    rename_bpm_arc(data_MQ_A,ia_K1L,fname,narc)

    for ii in range(0,nsec,1): 
        data_MQ_S=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.S{0}'.format(ii+1))]
        is_K1L=ref_tfs.columns.get_loc("K1L")
        rename_bpm_sec(data_MQ_S,is_K1L,fname,ii+1)
    
    fname.write('\n')
    
    fname.close()

#----------------------------------------------------------------------------------------------------------------

def svd_first(madx,i_start,seq,err,tol,beta):
    start=[]
    end=[]
    start.append('BPMN.A{0}.001'.format(i_start))
    if i_start%2 == 0:
        start.append('CORRN.S{0}.040'.format(i_start))
    else:
        start.append('CORRN.S{0}.052'.format(i_start))

    if i_start ==8:
        end.append('BPMN.S{0}.003'.format(1))
        end.append('CORRN.S{0}.001'.format(1))
    else:
        end.append('BPMN.S{0}.003'.format(i_start+1))
        end.append('CORRN.S{0}.001'.format(i_start+1))
        
    #end.append('CORRN.A{0}.369'.format(i_start))
 
    treading_svd(madx,i_start,start,end,seq,err,tol,beta)


def svd_ring(madx,seq,err,tol):
    start=[]
    end=[]
    for i in range(1,9,1):
        start.append('BPMN.A{0}.001'.format(i))
        if i%2 == 0:
            start.append('CORRN.S{0}.040'.format(i))
        else:
            start.append('CORRN.S{0}.052'.format(i))
        if i==8:    
            end.append('BPMN.S{0}.003'.format(1))
            end.append('CORRN.S{0}.001'.format(1))
        else:
            end.append('BPMN.S{0}.003'.format(i+1))
            end.append('CORRN.S{0}.001'.format(i+1))
        #end.append('CORRN.A{0}.369'.format(i))

    treading_svd_ring(madx,start,end,seq,err,tol)

#- - - - - -
'''    
def svd_last(madx,i_start,i_end,seq,err,tol,beta):
    start=[]
    end=[]
    start.append('IP{0}'.format(i_start))
    start.append('E.DIS_LEFT.A{0}'.format(i_start))
    end.append('FCC_HEBIP{0}_P_'.format(i_end))
    end.append('S.DIS_RIGHT.A{0}'.format(i_end))
    treading_svd(madx,i_start,start,end,seq,err,tol,beta)
'''
#- - - - - -

def treading_svd(madx,i_start,start,end,seq,err,tol,beta):
        madx.command.usemonitor(status='off')
        madx.command.usemonitor(status='on', range=start[0]+'/'+end[0])
        madx.command.usekick(status='off')
        madx.command.usekick(status='on', range=start[1]+'/'+end[1])

        #madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='x',corzero=0,clist='cx_fccee_heb_mic{0}.tab'.format(i_start))
        madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='x',corzero=0)

        madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='y',corzero=0)
        #madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='y',corzero=0,clist='cy_fccee_heb_mic{0}.tab'.format(i_start))

        twiss=madx.twiss(BETA0=beta, sequence=seq, file='FCCee_heb_modett_orbcor_arc{0}.tfs'.format(i_start))
        T=save_twiss(twiss)

        return

def treading_svd_ring(madx,start,end,seq,err,tol):
        madx.command.usemonitor(status='off')
        madx.command.usemonitor(status='on', range=start[0]+'/'+end[0])
        madx.command.usemonitor(status='on', range=start[2]+'/'+end[2])
        madx.command.usemonitor(status='on', range=start[4]+'/'+end[4])
        madx.command.usemonitor(status='on', range=start[6]+'/'+end[6])
        madx.command.usemonitor(status='on', range=start[8]+'/'+end[8])
        madx.command.usemonitor(status='on', range=start[10]+'/'+end[10])
        madx.command.usemonitor(status='on', range=start[12]+'/'+end[12])
        madx.command.usemonitor(status='on', range=start[14]+'/'+end[14])
        madx.command.usekick(status='off')
        madx.command.usekick(status='on', range=start[1]+'/'+end[1])
        madx.command.usekick(status='on', range=start[3]+'/'+end[3])
        madx.command.usekick(status='on', range=start[5]+'/'+end[5])
        madx.command.usekick(status='on', range=start[7]+'/'+end[7])
        madx.command.usekick(status='on', range=start[9]+'/'+end[9])
        madx.command.usekick(status='on', range=start[11]+'/'+end[11])
        madx.command.usekick(status='on', range=start[13]+'/'+end[13])
        madx.command.usekick(status='on', range=start[15]+'/'+end[15])

        madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='x',corzero=0,clist='cx_fccee_heb_mic_all.tab')
        madx.command.correct(flag='line',mode='svd',monerror=err,error=tol,plane='y',corzero=0,clist='cy_fccee_heb_mic_all.tab')

        return

#----------------------------------------------------------------------------------------------------------------
    
def rename_corr_arc(data,ind_k1,file_name,narc):
    ii=0
    for p in range(len(data)):
        jj=int(p%(len(data)/narc)+1)
        if jj==1:
            ii+=1
            
        if data.iloc[p,ind_k1]>0:
            if jj<10:
                file_name.write("corrn.a{0}.00{1} : hcorr,kick:=k.mcbh.a{0}.00{1};\n".format(ii,jj))
            elif jj>=10 and jj<100:
                file_name.write("corrn.a{0}.0{1} : hcorr,kick:=k.mcbh.a{0}.0{1};\n".format(ii,jj))
            elif jj>=100 and jj<370:
                file_name.write("corrn.a{0}.{1} : hcorr,kick:=k.mcbh.a{0}.{1};\n".format(ii,jj))
                
        elif data.iloc[p,ind_k1]<0:
            if jj<10:
                file_name.write("corrn.a{0}.00{1} : vcorr,kick:=k.mcbv.a{0}.00{1};\n".format(ii,jj))
            elif jj>=10 and jj<100:
                file_name.write("corrn.a{0}.0{1} : vcorr,kick:=k.mcbv.a{0}.0{1};\n".format(ii,jj))
            elif jj>=100 and jj<370:
                file_name.write("corrn.a{0}.{1} : vcorr,kick:=k.mcbv.a{0}.{1};\n".format(ii,jj))

#- - - - - -
def rename_corr_sec(data,ind_k1,file_name,nsec):
    for p in range(len(data)):
        jj=int(p%len(data)+1)

        if data.iloc[p,ind_k1]>0:
            if jj<10:
                file_name.write("corrn.s{0}.00{1} : hcorr,kick:=k.mcbh.s{0}.00{1};\n".format(nsec,jj))
            elif jj>=10 and jj<100:
                file_name.write("corrn.s{0}.0{1} : hcorr,kick:=k.mcbh.s{0}.0{1};\n".format(nsec,jj))

        elif data.iloc[p,ind_k1]<0:
            if jj<10:
                file_name.write("corrn.s{0}.00{1} : vcorr,kick:=k.mcbv.s{0}.00{1};\n".format(nsec,jj))
            elif jj>=10 and jj<100:
                file_name.write("corrn.s{0}.0{1} : vcorr,kick:=k.mcbv.s{0}.0{1};\n".format(nsec,jj))
        
#- - - - - -
                
def replace_corr(madx,file_ref,narc,nsec):

    fnameref=file_ref
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    ref_tfs=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

    data_MQ_A=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A')]
    
    ii=0
    for p in range(len(data_MQ_A)):
        jj=int(p%(len(data_MQ_A)/narc)+1)
        if jj==1:
            ii+=1

        if jj<10:
            madx.replace(element='corr.a{0}.00{1}'.format(ii,jj), by='corrn.a{0}.00{1}'.format(ii,jj))
        elif jj>=10 and jj<100:
            madx.replace(element='corr.a{0}.0{1}'.format(ii,jj), by='corrn.a{0}.0{1}'.format(ii,jj))
        elif jj>=100 and jj<370:
            madx.replace(element='corr.a{0}.{1}'.format(ii,jj), by='corrn.a{0}.{1}'.format(ii,jj))

    for ii in range(0,nsec,1): 
        data_MQ_S=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.S{0}'.format(ii+1))]
    
        #ii=0
        for p in range(len(data_MQ_S)):
            jj=int(p%len(data_MQ_S)+1)
            if jj==1:
                ii+=1

            if jj<10:
                madx.replace(element='corr.s{0}.00{1}'.format(ii,jj), by='corrn.s{0}.00{1}'.format(ii,jj))
            elif jj>=10 and jj<100:
                madx.replace(element='corr.s{0}.0{1}'.format(ii,jj), by='corrn.s{0}.0{1}'.format(ii,jj))
            elif jj>=100 and jj<370:
                madx.replace(element='corr.s{0}.{1}'.format(ii,jj), by='corrn.s{0}.{1}'.format(ii,jj))
                
#- - - - - -

def macro_corr(file_ref,file_out,narc,nsec):

    fnameref=file_ref
    head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
    ref_tfs=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)

    if os.path.exists(file_out):
        os.remove(file_out)
    else:
        print('Not possible: the file does not exist')

    fname=open(file_out, "w+")
    
    #fname.write("HVCORR : macro ={\n")

    fname.write("hcorr : hkicker, l:=l_corr;\n")
    fname.write("vcorr : vkicker, l:=l_corr;\n")

    data_MQ_A=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.A')]
    ia_K1L=ref_tfs.columns.get_loc("K1L")
    rename_corr_arc(data_MQ_A,ia_K1L,fname,narc)
    
    for ii in range(0,nsec,1): 
        data_MQ_S=ref_tfs.loc[ref_tfs["NAME"].str.contains('MQ.S{0}'.format(ii+1))]
        is_K1L=ref_tfs.columns.get_loc("K1L")
        rename_corr_sec(data_MQ_S,is_K1L,fname,ii+1)
        
    fname.write('\n')
    
    fname.close()

#----------------------------------------------------------------------------------------------------------------
    
