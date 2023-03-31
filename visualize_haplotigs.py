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
len_cutoff = 5000

cmap = cmr.neon
hap = argv[1]
hap_cov = []
lines = []
line_colors = []
cov_cutoff = 1.5
widths = []
max_br = 0
cont = False
for line in open(hap,'r'):
    if line[0] == '>':
        name = line.split()[0][1:]
        ar = p.findall(line)
        br = snp_p.findall(line)
        if len(ar) == 0:
            cont = True
            print(line)
            continue
        curr_cov = float(ar[0])
        start = int(br[0][0])
        end = int(br[0][1])
        if end - start < len_cutoff:
            cont = True
            continue
        baserange = (start,end)
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
            lines.append([(baserange[0], (curr_cov+1)), (baserange[1], (curr_cov+1))])
            line_colors.append(total_alt / al_c)
            widths.append(4)

lc = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors), cmap = cmap, alpha = 1.0)
print(lc)

fig, ax = plt.subplots()
fig.set_size_inches(8 * cm, 6 * cm)
ax.add_collection(lc)
ax.set_xlim(0, max_br)
ax.set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel("Base position")
plt.ylabel("Average coverage of vartig")
t = hap.split('-')[0]
plt.title(t)
#fig.colorbar(lc, label = "Alternate allele ratio")

plt.tight_layout()
plt.savefig(f"vartig-plot_{hap}.png", dpi = 300)
plt.show()

