#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import pathlib


# function to add value labelsdef addlabels(x,y):
def add_value_label(x_list,y_list):
    for i in range(1, len(x_list)+1):
        plt.text(i,y_list[i-1]/2,y_list[i-1], ha="center")


l_err_mq=[20,40,60,80,90,100]
eseed=100
l_passed=[]

for i in range(len(l_err_mq)):
    cnt=0
    for j in range(eseed):
        fname='/home/td271008/work/cpymadx/Orbit_cpymad/mq_offset_{0}_IP5_2it_{1}seeeds_corrhplus/test_seed{2}/FCCee_heb_modett_orbcor_all_sextuon_it2_seed{2}.tfs'.format(l_err_mq[i],eseed,j+1)

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

fig.savefig('/home/td271008/work/cpymadx/Orbit_cpymad/comp_mq_offset.pdf')
#plt.show()
