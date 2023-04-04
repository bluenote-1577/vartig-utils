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



threshold = [0.00,1.00]
mult = 1.0
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
cm = 1/2.54  # centimeters in inches
chunk = 200

cmap = cmr.tropical
haps = argv[1:]
hap_cov = []
lines = []
line_colors = []
cov_cutoff = 1.5
widths = []
max_br = 0
cont = False
fig, ax = plt.subplots()
for hap in haps:
    running_c = 0
    total_alt = 0
    num_zero = 0
    last_pos = 0

    for line in open(hap,'r'):
        if line[0] == '>':
            continue
        else:
            spl = line.split()
            if running_c == 0:
                last_pos = int(spl[0].split(':')[1])
            vals = spl[2].split('|')
            for thing in vals:
                tup  = thing.split(':')
                if tup[0] == "NA":
                    running_c += 1
                    continue
                al = int(tup[0])
                count = int(tup[1])
                if al == 0:
                    num_zero += count
                else:
                    total_alt += count

            running_c += 1
            if running_c == chunk:
                running_c = 0
                cur_pos = int(spl[0].split(':')[1])
                cur_cov = num_zero + total_alt
                lines.append([(last_pos, cur_cov / chunk), (cur_pos, cur_cov/chunk )])

                print(last_pos,cur_pos,hap )
                last_pos = cur_pos

                hap_cov.append(cur_cov/chunk)
                line_colors.append(total_alt / cur_cov)
                widths.append(4)
                num_zero = 0
                total_alt = 0
                if last_pos > max_br:
                    max_br = last_pos


lc = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors), cmap = cmap, alpha = 1.0)

ax.add_collection(lc)
ax.set_xlim(0, max_br)
ax.set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)

plt.show()

