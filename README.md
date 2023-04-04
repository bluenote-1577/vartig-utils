# vartig-utils - utilities for processing and aligning vartigs between samples


## TODO INTERNAL - JIM

This is the repo I used for generating strain tracking information and visualize vartigs. The "vartig" files below mean the glopp-out-dir/contig/haplotigs.fa file. 

## Introduction

This is a companion repository for processing outputs from our software TODO. We currently support mapping of vartigs between samples
(e.g. to track strains longitudinally) or computing statistics between vartigs. 

This repo also contains scripts for visualizing vartigs.

### Installation

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)

For the python scripts:

1. numpy
2. scipy
3. matplotlib.
4. [cmasher](https://cmasher.readthedocs.io/user/introduction.html#how-to-install). This could be skipped, but you'll have to change the colormap used in the scripts.

Building takes a few minutes (depending on # of cores).

```sh
git clone https://github.com/bluenote-1577/vartig-utils
cd vartig-utils

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
vtig -h

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
#./target/release/vtig -h
```

## Quick start

```sh
# compare two vartigs and display some statistics 
vtig dist vartig1 vartig2

# find best matching vartigs between two samples
vtig map vartig1 vartig2 -m minimum_match_length

# visualize vartigs
python visualize_vartigs.py vartig1

# track vartigs between samples. REQUIRES `vtig` command to be in path. 
python track_vartigs.py CONTIG_NAME vartig1 vartig2 vartig3

```
All above scripts/commands work **only if** vartigs are generated from **the exact same vcf and reference**. 

You SHOULD NOT use this for:

1. mapping vartigs from different contigs (even from the same MAG, for example) 
2. mapping vartigs from the same contigs, but a different VCF (you must use the *same* vcf to generate both vartigs)


## Output

The output for `vtig map` looks like:
```
name1   name2   identity        num_same_allele num_diff_allele cov1    cov2    snp_range1      snp_range2      base_range1     base_range2
sample1/NZ_AP024085.1_HAP15   sample2/NZ_AP024085.1_HAP13 1.000       7       0       224.28      224.28      17-23   17-23   9960-10786      9960-10786
sample1/NZ_AP024085.1_HAP16   sample2/NZ_AP024085.1_HAP14 1.000       7       0       39.71      39.71     17-23   17-23   9960-10786      9960-10786

```
- name1: the vartig name of the first sample
- name2: the vartig name of the second sample
- identity: num_same_allele / (num_same_allele + num_diff_allele)
- num_same_allele/num_diff_allele: number of matching and differing alleles between the two vartigs
- cov1/cov2: the average coverage over all alleles. 
- snp_range1/2: the range of SNPs covered by the vartig
- base_range1/2: the range of bases on the contig covered by the vartig

The output for `python visualize_vartigs.py` (for a real, nanopore phasing of *Brevefilum fermentans*) looks like this:

![Visualization example](https://github.com/bluenote-1577/vartig-utils/blob/main/visualize-vartig-example.png)

Each bar represents a vartig in the file. 

- The color in the top plot represents the fraction of alternate alleles in the vartig. For example, if the vartig has all alternate SNPs, then it will have 1.0 alternate allele fraction.
- The color in the bottom plot represents the HAPQ, which is a heuristic measure of phasing ambiguity between 0-60, analogous to MAPQ. Notice the small blocks have low HAPQ, most likely indicating they're spurious. 


## Citation

TODO
