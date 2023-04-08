import matplotlib.pyplot as plt
import cmasher as cmr
def normal (x):
    return np.log2(x+np.array(1))


from scipy import stats

import numpy as np
from matplotlib import collections  as mc

import os
import re
from sys import argv
import subprocess
import shlex
p = re.compile('COV:(\d*\.?\d+)')
snp_p = re.compile('BASERANGE:(\d+)-(\d+)')
hapq_p = re.compile('HAPQ:(\d+)')
if len(argv) > 2:
    hapq_cut = int(argv[2])

else:
    hapq_cut = 0


threshold = [0.00,1.00]
mult = 1.0
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
cm = 1/2.54  # centimeters in inches
len_cutoff = 5000

cmap = cmr.neon
#cmap = plt.cm.viridis_r
cmap_hq = plt.cm.binary
hap = argv[1]
hap_cov = []
lines = []
line_colors = []
line_colors_hapq = []
cov_cutoff = 1.5
widths = []
max_br = 0
cont = False
hapq_fail = False
for line in open(hap,'r'):
    if line[0] == '>':
        name = line.split()[0][1:]
        ar = p.findall(line)
        br = snp_p.findall(line)
        hapr = hapq_p.findall(line)
        if len(ar) == 0:
            cont = True
            print(line)
            continue
        curr_cov = float(ar[0])
        start = int(br[0][0])
        end = int(br[0][1])
        hapq = int(hapr[0])
        if hapq > 60:
            hapq = 60
        print(hapq)
        baserange = (start,end)
        if end - start < len_cutoff:
            cont = True
            continue
        if hapq < hapq_cut:
            hapq_fail = True
        if end > max_br:
            max_br = end
    else:
        if cont:
            cont = False
            continue
        total_alt = 0
        al_c = 0
        num_zero = 0
        for x in line:
            if x == '0' or x == '1' or x == '2':
                al_c += 1
            if x == '1' or x == '2':
                total_alt += 1 
            else:
                num_zero += 1
        if curr_cov > cov_cutoff:
            if al_c == 0:
                print(line, curr_cov)
                continue
            hap_cov.append(curr_cov)
            lines.append([(baserange[0], (curr_cov+1)), (baserange[1] + 5000, (curr_cov+1))])
            line_colors.append(total_alt / al_c)
            line_colors_hapq.append(hapq)
            if hapq_fail:
                hapq_fail = False
                widths.append(1)
            else:
                widths.append(4)
        else:
            hapq_fail = False

lc = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors), cmap = cmap, alpha = 1.0)
lc_hq = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors_hapq), cmap = cmap_hq, alpha = 1.0)
lc_hq.set_clim(vmin=0, vmax=60)
lc.set_clim(vmin=0, vmax=1)
print(lc)

fig, ax = plt.subplots(2)
ax[0].add_collection(lc)
ax[0].set_xlim(0, max_br)
ax[0].set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)
fig.colorbar(lc, ax=ax[0], orientation='vertical')
fig.colorbar(lc_hq, ax=ax[1], orientation='vertical')



ax[1].add_collection(lc_hq)
ax[1].set_xlim(0, max_br)
ax[1].set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)
ax[0].set_title("Colored by alternate allele ratio")
ax[1].set_title("Colored by HAPQ")
plt.show()

