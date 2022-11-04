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


#l_err_mq=[60,80,90,100,120,150,200] #with mq offsets, mb err field & roll
l_err_mq=[60,100,150,200] #with mq offsets, mb err field & roll + tm
eseed=100
l_passed=[]

Q1=415.225
Q2=416.29

#if 1: compares without tune match i.e. if twiss in ring with sextuon passes
#if 2: check if the match of the tune is well done
FLAG=1


if FLAG==1:
    for i in range(len(l_err_mq)):
        cnt=0
        for j in range(eseed):

            fname=path+'mb_b1r_dpsi_mq_dxdy_{0}_{1}seeds_corrhplus/test_seed{2}/FCCee_heb_modett_orbcor_all_ring_sextuon_it1_seed{2}.tfs'.format(l_err_mq[i],eseed,j+1)

            if os.path.exists(fname):
                cnt+=1
            else:
                continue

        l_passed.append((cnt*eseed)/100)


    x = np.arange(len(l_err_mq))
    width=0.35
    fig, ax = plt.subplots()

    plt.xlabel("Offset MQ ($\mu$m)")
    plt.ylabel("% of successful seeds")
    ax.set_xticks(x)
    ax.set_xticklabels(l_err_mq)

    pps = ax.bar(x, l_passed, width)
    for p in pps:
       height = p.get_height()
       ax.annotate('{}'.format(height),
          xy=(p.get_x() + p.get_width() / 2, height),
          xytext=(0, 3), # 3 points vertical offset
          textcoords="offset points",
          ha='center', va='bottom')

    fig.savefig(path+'comp_mb_fielderr_roll_mq_offset.pdf')
    #plt.show()




if FLAG==2:
    for i in range(len(l_err_mq)):
        cnt=0
        for j in range(eseed):

            fname=path+'mb_b1r_dpsi_mq_dxdy_{0}_tm_{1}seeds_corrhplus/test_seed{2}/FCCee_heb_modett_tune_match_seed{2}.tfs'.format(l_err_mq[i],eseed,j+1)

            if os.path.exists(fname):
                head_opt=pd.read_csv(fname, header=50, sep='\s+', nrows=0).columns[1:]
                headers=pd.read_csv(fname, sep='\s+', names=head_opt, on_bad_lines='skip', low_memory=False)
                Q_x=float(headers.iat[22,3]) #tune in plane x
                Q_y=float(headers.iat[23,3]) #tune in plane y

                if Q_x==Q1 and Q_y==Q2:
                    cnt+=1
            else:
                continue

        l_passed.append((cnt*eseed)/100)


    x = np.arange(len(l_err_mq))
    width=0.35
    fig, ax = plt.subplots()

    plt.xlabel("Offset MQ ($\mu$m)")
    plt.ylabel("% of successful seeds")
    ax.set_xticks(x)
    ax.set_xticklabels(l_err_mq)

    pps = ax.bar(x, l_passed, width)
    for p in pps:
       height = p.get_height()
       ax.annotate('{}'.format(height),
          xy=(p.get_x() + p.get_width() / 2, height),
          xytext=(0, 3), # 3 points vertical offset
          textcoords="offset points",
          ha='center', va='bottom')

    fig.savefig(path+'comp_tm_mb_fielderr_roll_mq_offset.pdf')
