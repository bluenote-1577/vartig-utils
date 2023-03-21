import matplotlib.pyplot as plt
import os
import re
from sys import argv
import subprocess
import shlex

haps = argv[1:]

hap_covs = []
len_cutoff = 2
vtig_cutoff = 5

p = re.compile('COV:(\d*\.?\d+)')
for hap in haps:
    hap_covs.append([])
    curr_cov = 0
    for line in open(hap,'r'):
        if line[0] == '>':
            ar = p.findall(line)
            curr_cov = float(ar[0])
        else:
            al_c = 0
            for x in line:
                if x == '0' or x == '1' or x == '2':
                    al_c += 1
            if al_c > vtig_cutoff:
                hap_covs[-1].append(curr_cov)


haps = [shlex.quote(x) for x in haps]
for i in range(len(haps)-1):

    x = []
    y = []
    count_miss = 0
    count = 0

    os.system(f"vtig map {haps[i]} {haps[i+1]} -m {len_cutoff} > x")
    for line in open("x",'r'):
        if line[0] == "#":
            continue
        spl = line.split()
        x.append(float(spl[5]))
        y.append(float(spl[6]))

        if float(spl[2]) > 0.95 and float(spl[4]) <= 1:
            tt = 1
            plt.plot([i,i+1],[x[-1],y[-1]], 'b', alpha = min(float(spl[3]) / 1000,1))
            plt.plot([i,i+1],[x[-1],y[-1]], 'b', alpha = float(0.01))
        else:
            count_miss += 1
        count += 1



    print(count_miss, count)
    ones = [i for x in range(len(hap_covs[i]))]
    twos  = [i+1 for x in range(len(hap_covs[i+1]))]
    plt.plot(ones,hap_covs[i],'o', ms = 4, c = 'b', alpha = 0.05)
    plt.plot(twos,hap_covs[i+1],'o', ms = 4, c = 'b', alpha = 0.05)

plt.ylabel("coverage")
plt.title("haplotig correspondence from sample1 to sample 2 (each line is haplotig)")
plt.xlabel("sample")

plt.show()

