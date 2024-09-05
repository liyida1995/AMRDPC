# AMRDPC: An Efficient Lossy Compression Framework for Density Partitioning in AMR Applications.

## Overview
The main idea of AMRDPC is to preprocess AMR data based on density partitioning before compression, and further improve AMR applications' storage efficiency.

## Instruction for AMRDPC
1. Get the AMR applications and lossy compressor SZ.
2.

## Evaluation Platform
We conduct our experiments on a Linux server with OS kernel of Linux 5.15, CPU of 12th Gen Intel(R) Core(TM) i9-12900K, main memory of 32 GB DDR5 RAM, and storage device of 1TB M.2 Gen4 NVMe SSD.

## AMR Applications
We evaluate the effectiveness of AMRDPC with seven real AMR application datasets in AMReX. These applications can be downloaded from the following addresses:
1. Grid Z2: https://drive.google.com/file/d/118ataZqVlCFq5D04MhtTTOBBvVO6DR6B/view?usp=sharing
2. Grid Z3: https://drive.google.com/file/d/1PHMUV6c9K5V9AN-TF_CYf9ROR86oZYSm/view?usp=sharing
3. Grid Z5: https://drive.google.com/file/d/1hEjnw3llHj2P9z3WyrWSQAWqPAQ-LKxy/view?usp=sharing
4. Grid Z10: https://drive.google.com/file/d/1xWwf2_YFDjod0mUIRl1iffj95ZemR4Q6/view?usp=sharing
5. Run2_T2: https://drive.google.com/file/d/1xZnts_QRAN5Av3NbBzLzt0FybVafV7JC/view?usp=sharing
6. Run2_T3: https://drive.google.com/file/d/1vy3Bry1RmAdXqMNXhwXFWZOzF6bCa7hM/view?usp=sharing
7. Run2_T4: https://drive.google.com/file/d/1CeI_nxhU5_rrT2e78Jsy025t5lOmaQAd/view?usp=sharing

## Data Compression
We utilize the lossy compressors SZ to evaluate how much AMRDPC can improve the compression ratio. In particular, we use SZ with the default mode and absolute error bound. SZ can be downloaded from the following address:
1. SZ: https://github.com/szcompressor/SZ.

## Comparison Baseline
1. The primitive 3D baseline: Different AMR levels are unified to the same resolution for 3D compression;
2. The state-of-the-art 3D compressor TAC, which involves multiple density-based strategies to process AMR data.

## Testing Examples
Examples can be found in the AMRDPC/example.

You can use the executable 'sh test.sh' command to do the compression.


## Overall Results
Experimental results show that compared with the state-of-the-art method, AMRDPC improves the data compression ratio by up to 5.73X while ensuring high data fidelity. The data compression speed is increased by up to 18.83% without reducing the data quality.
