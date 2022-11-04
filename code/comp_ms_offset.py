#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import pathlib

path="/feynman/work/dacm/leda/td271008/Orbit_cpymad/"


# function to add value labelsdef addlabels(x,y):
def add_value_label(x_list,y_list):
    for i in range(1, len(x_list)+1):
        plt.text(i,y_list[i-1]/2,y_list[i-1], ha="center")



l_err_ms=[60,80,100,120,150,200] #offset ms
eseed=100
l_passed=[]

#if 1: compares without tune match i.e. if twiss in ring with sextuon passes
#if 2: check if the match of the tune is well done (without svd)
#if 3: check if the match of the tune is well done (with svd)
FLAG=3


if FLAG==1:
    print('FLAG1')
    for i in range(len(l_err_ms)):
        cnt=0
        for j in range(eseed):

            fname=path+'ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_{1}seeds/test_seed{2}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{2}.tfs'.format(l_err_ms[i],eseed,j+1)

            if os.path.exists(fname):
                cnt+=1
            else:
                continue

        l_passed.append((cnt*eseed)/100)


    x = np.arange(len(l_err_ms))
    width=0.35
    fig, ax = plt.subplots()

    plt.xlabel("Offset MS ($\mu$m)")
    plt.ylabel("% of successful seeds")
    ax.set_xticks(x)
    ax.set_xticklabels(l_err_ms)

    pps = ax.bar(x, l_passed, width)
    for p in pps:
       height = p.get_height()
       ax.annotate('{}'.format(height),
          xy=(p.get_x() + p.get_width() / 2, height),
          xytext=(0, 3), # 3 points vertical offset
          textcoords="offset points",
          ha='center', va='bottom')

    fig.savefig(path+'comp_ms_offset_bpm_offset_mb_fielderr_roll_mq_offset.pdf')
    #plt.show()




if FLAG==2:
    print('FLAG2')
    for i in range(len(l_err_ms)):
        cnt=0
        for j in range(eseed):

            file1=path+'ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_{1}seeds/test_seed{2}/FCCee_heb_modett_seed{2}.tfs'.format(l_err_ms[i],eseed,j+1)
            head_opt=pd.read_csv(file1, header=50, sep='\s+', nrows=0).columns[1:]
            headers=pd.read_csv(file1, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
            Q1=float(headers.iat[22,3]) #tune in plane x
            Q2=float(headers.iat[23,3]) #tune in plane y

            fname=path+'ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_{1}seeds/test_seed{2}/FCCee_heb_modett_tune_match_it0_seed{2}.tfs'.format(l_err_ms[i],eseed,j+1)

            if os.path.exists(fname):
                head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
                headers=pd.read_csv(fname, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
                Q_x=float(headers.iat[22,3]) #tune in plane x
                Q_y=float(headers.iat[23,3]) #tune in plane y

                if round(Q_x,3)==round(Q1,3) and round(Q_y,2)==round(Q2,2):
                    cnt+=1
            else:
                continue

        l_passed.append((cnt*eseed)/100)


    x = np.arange(len(l_err_ms))
    width=0.35
    fig, ax = plt.subplots()

    plt.xlabel("Offset MS ($\mu$m)")
    plt.ylabel("% of successful seeds")
    ax.set_xticks(x)
    ax.set_xticklabels(l_err_ms)

    pps = ax.bar(x, l_passed, width)
    for p in pps:
       height = p.get_height()
       ax.annotate('{}'.format(height),
          xy=(p.get_x() + p.get_width() / 2, height),
          xytext=(0, 3), # 3 points vertical offset
          textcoords="offset points",
          ha='center', va='bottom')

    fig.savefig(path+'comp_tm_it0_ms_offset_bpm_offset_mb_fielderr_roll_mq_offset.pdf')



if FLAG==3:
    print('FLAG3')
    for i in range(len(l_err_ms)):
        cnt=0
        for j in range(eseed):

            file1=path+'ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_{1}seeds/test_seed{2}/FCCee_heb_modett_seed{2}.tfs'.format(l_err_ms[i],eseed,j+1)
            head_opt=pd.read_csv(file1, header=50, sep='\s+', nrows=0).columns[1:]
            headers=pd.read_csv(file1, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
            Q1=float(headers.iat[22,3]) #tune in plane x
            Q2=float(headers.iat[23,3]) #tune in plane y

            fname=path+'ms_dxdy_{0}_bpm_dxdy_150_mb_b1r_dpsi_mq_dxdy_150_tm_svd_{1}seeds/test_seed{2}/FCCee_heb_modett_tune_match_it1_seed{2}.tfs'.format(l_err_ms[i],eseed,j+1)

            if os.path.exists(fname):
                head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
                headers=pd.read_csv(fname, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
                Q_x=float(headers.iat[22,3]) #tune in plane x
                Q_y=float(headers.iat[23,3]) #tune in plane y

                if round(Q_x,3)==round(Q1,3) and round(Q_y,2)==round(Q2,2):
                    cnt+=1
            else:
                continue

        l_passed.append((cnt*eseed)/100)


    x = np.arange(len(l_err_ms))
    width=0.35
    fig, ax = plt.subplots()

    plt.xlabel("Offset MS ($\mu$m)")
    plt.ylabel("% of successful seeds")
    ax.set_xticks(x)
    ax.set_xticklabels(l_err_ms)

    pps = ax.bar(x, l_passed, width)
    for p in pps:
       height = p.get_height()
       ax.annotate('{}'.format(height),
          xy=(p.get_x() + p.get_width() / 2, height),
          xytext=(0, 3), # 3 points vertical offset
          textcoords="offset points",
          ha='center', va='bottom')

    fig.savefig(path+'comp_tm_it1_ms_offset_bpm_offset_mb_fielderr_roll_mq_offset.pdf')
