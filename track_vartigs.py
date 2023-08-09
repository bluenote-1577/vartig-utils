import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#figure order
#$CONTIG haplotigs/glopp-SRR13712920-ab-filter_{$CONTIG}_haplotigs.fa   haplotigs/glopp-SRR13712919-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712896-ab-filter_{$CONTIG}_haplotigs.fa  haplotigs/glopp-SRR13712866-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712899-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712885-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712854-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712844-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712833-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712822-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712918-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712908-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712907-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712905-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712904-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712903-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712906-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712902-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712901-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712900-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712897-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712874-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712873-ab-filter_{$CONTIG}_haplotigs.fa haplotigs/glopp-SRR13712898-ab-filter_{$CONTIG}_haplotigs.fa


import cmasher as cmr
from scipy import stats

import numpy as np
from matplotlib import collections  as mc

def normal (x, normal_factor):
    if normal_factor == False:
        #return 2 * np.sqrt(x + np.array(3/8))
        return np.log2(x+np.array(1))
        #return np.sqrt(x)
    else:
        return x/normal_factor


threshold = [0.25,0.75]
mult = 1.0
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
cm = 1/2.54  # centimeters in inches

#cmap=plt.cm.PiYG
#cmap = plt.cm.viridis_r

#cmap = cmr.redshift
#cmap = cmr.watermelon
cmap = cmr.neon
#cmap = cmr.tropical
#cmap = cmr.pride
#cmap = cmr.bubblegum_r
#cmap = cmr.iceburn
#cmap=plt.cm.cividis


import os
import re
from sys import argv
import subprocess
import shlex
mean_normal = False
line_offset = 0.04

contig = argv[1]
haps = argv[2:]
kde_weights_none = [None for _ in haps]
avg_alt = 0
avg_al_c = 0

hap_covs = []
kde_weights = []
len_cutoff = 2
cov_cutoff = 1.5
dates = [0,3,4,5,41,48,49,301,339,346,354,360,391,446,452,473,476,516,542,553,614,622,627,636]
vartig_to_allele_frac = dict()
print(len(argv))

p = re.compile('COV:(\d*\.?\d+)')
snp_p = re.compile('SNPRANGE:(\d+)-(\d+)')
fig,ax = plt.subplots(2)
fig.set_size_inches(17 * cm, 10 * cm)
for hap in haps:
    print(hap)
    hap_covs.append([])
    kde_weights.append([])
    curr_cov = 0
    for line in open(hap,'r'):
        if line[0] == '>':
            name = line.split()[0][1:]
            vartig_to_allele_frac[name] = None
            ar = p.findall(line)
            curr_cov = float(ar[0])
        else:
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
            if al_c > len_cutoff and al_c > cov_cutoff:
                hap_covs[-1].append(curr_cov)
                kde_weights[-1].append(al_c)
                vartig_to_allele_frac[name] = [al_c, total_alt]
                avg_alt += total_alt
                avg_al_c += al_c


avg_minor = avg_alt / avg_al_c
#print(avg_minor)

haps = [shlex.quote(x) for x in haps]

hcs = [[] for _ in haps]

for i in range(len(haps)-1):
    lines = []
    widths = []
    x = []
    y = []
    count_miss = 0
    count = 0
    hc_1 = hcs[i]
    hc_2 = hcs[i+1]
    line_colors = []

    os.system(f"vtig map {haps[i]} {haps[i+1]} -m {len_cutoff} > x")
    os.system(f"vtig dist {haps[i]} {haps[i+1]} -m {len_cutoff}")
    for line in open("x",'r'):
        #first line
        if 'name1' in line:
            continue

        if line[0] == "#":
            continue
        spl = line.split()
        if float(spl[5]) > 500 or float(spl[6]) > 500:
            continue
        
        if float(spl[2]) > 0.95 and float(spl[4]) <= 1 and float(spl[5]) > cov_cutoff:
            x.append(float(spl[5]))
            y.append(float(spl[6]))
            hc_1.append(x[-1])
            hc_2.append(y[-1])
            lines.append([(i,x[-1]),(i+1,y[-1])])
            vals1 =  vartig_to_allele_frac[spl[0]]
            vals2 =  vartig_to_allele_frac[spl[1]]
            if vals1 is None or vals2 is None:
                continue
            val = (vals1[1]+ vals2[1])/(vals1[0] + vals2[0])
            #if val > avg_minor:
            #    val = 1
            #else: 
            #    val = 0
            #adj = 0.5 - avg_minor
            adj = 0
            val = (val + adj) * mult;
            line_colors.append(np.max([np.min([val, threshold[1]]), threshold[0]]))
            tt = 1
            #widths.append(float(spl[3])/500)
            widths.append(1/50)
            #plt.plot([i,i+1],[x[-1],y[-1]], 'b', alpha = min(float(spl[3]) / 1000,1))
            #plt.plot([i,i+1],[x[-1],y[-1]], 'b', alpha = float(0.01))
        else:
            count_miss += 1
        count += 1


    line_colors = np.array(line_colors)
    mean_x = np.mean(hap_covs[i]) 
    mean_y = np.mean(hap_covs[i+1])
    #mean_x = 1
    #mean_y = 1
    #lines_normal = [[(i, x[j]/mean_x),(i+1,y[j]/mean_y)] for j in range(len(x))]
    lines_normal = [[(i+line_offset, normal(x[j], mean_normal)),(i+1 - line_offset,normal(y[j], mean_normal))] for j in range(len(x))]
    lines = [[(i+line_offset, x[j]),(i+1-line_offset,y[j])] for j in range(len(x))]
    lc = mc.LineCollection(lines, linewidths=widths, array = line_colors, cmap = cmap, alpha = 0.5)
    lc_normal = mc.LineCollection(lines_normal, linewidths=widths, array = line_colors, cmap = cmap, alpha = 0.5)

    #lc = mc.LineCollection(lines, linewidths=widths, array = line_colors)
    #lc_normal = mc.LineCollection(lines_normal, linewidths=widths, array = line_colors)

    ax[0].add_collection(lc)
    ax[1].add_collection(lc_normal)
    #print(count_miss, count)
        #ax[0].plot(ones,np.array(hap_covs[i]),'o', ms = ms, c = 'b', alpha = al)
    #ax[0].plot(twos,np.array(hap_covs[i+1]),'o', ms = ms, c = 'b', alpha = al)
    #ax[1].plot(ones,np.array(normal(hap_covs[i], mean_normal)),'o', ms = ms, c = 'b', alpha = al)
    #ax[1].plot(twos,np.array(normal(hap_covs[i+1], mean_normal)),'o', ms = ms, c = 'b', alpha = al)
    

for i in range(len(hap_covs) - 1):

    weights = None
    ones = [i for x in range(len(hap_covs[i]))]
    twos  = [i+1 for x in range(len(hap_covs[i+1]))]
    al = 0.0020
    ms = 1
    bw = None

    #print(len(kde_weights[i]))
    #print(len(hap_covs[i]))
    norm_vals = normal(hap_covs[i], mean_normal)
    kde_norm = stats.gaussian_kde(norm_vals, bw_method = bw, weights = kde_weights[i])
    kde_norm.bw_method = bw
    c_norm = kde_norm(norm_vals)

    #plt.hist(normal(hap_covs[i],False),bins=200, density = True)
    #plt.plot(norm_vals, c_norm, 'o')
    #plt.show()
    #exit()


    kde = stats.gaussian_kde(hap_covs[i], weights = kde_weights[i])
    kde.bw_method = bw
    c = kde(hap_covs[i])

    print(c_norm/np.max(c_norm))

    ax[1].scatter(ones, norm_vals, c = c_norm, cmap = "binary", s = ms, alpha = 1)
    ax[0].scatter(ones, hap_covs[i], c = c, cmap = "binary", s = ms, alpha = 1)

    if i == len(hap_covs) - 2:
        norm_vals = normal(hap_covs[i+1], mean_normal)
        kde_norm = stats.gaussian_kde(norm_vals, weights = kde_weights[i+1])
        kde_norm.bw_method = bw
        c_norm = kde_norm(norm_vals)

        kde = stats.gaussian_kde(hap_covs[i+1], weights = kde_weights[i+1])
        kde.bw_method = bw
        c = kde(hap_covs[i+1])

        ax[1].scatter(twos, norm_vals, c = c_norm, cmap = "binary", s = ms, alpha = 1)
        ax[0].scatter(twos, hap_covs[i+1], c = c, cmap = "binary", s = ms, alpha = 1)


if len(hcs) == len(dates) and ("NZ_AP024085.1" in contig or "NC_021016.1" in contig):
    ax[0].set_xticks([])
    plt.xticks(range(len(hcs)),dates, rotation='vertical')


ymax = np.max([np.percentile(hap_covs[i], 99) for i in range(len(hap_covs))])
ymax_norm = np.max([np.percentile(normal(hap_covs[i],mean_normal), 99) for i in range(len(hap_covs))])
ax[0].set_ylim([0,ymax])
ax[1].set_ylim([0,ymax_norm])
ax[0].set_ylabel("Coverage")
ax[1].set_ylabel("Log2(coverage+1)")
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)

if "NZ_AP024085.1" in contig:
    title = f"Faecalibacillus intestinalis"
elif "NC_021016.1":
    title = f"Anaerostipes hadrus"
else:
    title = f"{contig}"

ax[0].set_title(title)
plt.xlabel("Days since first sample")
#ax[2].axis('off')
fig.colorbar(lc, ax = ax[0], label = "Alternate allele ratio")
fig.colorbar(lc, ax = ax[1], label = "Alternate allele ratio")
#fig.colorbar(lc,  ax=ax[1])

plt.tight_layout()
plt.savefig(f"24_track_{contig}.png", dpi = 300)
plt.show()

