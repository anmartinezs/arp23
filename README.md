# arp23
Python scripts to statistically describe Arp2/3 complex organization on travelling Actin Waves

## Contents and usage
This project contains the scripts used to generate some of the results published in:

> Jasnin et al. (2019) The Architecture of Traveling Actin Waves Revealed by Cryo-Electron Tomography. Accepted in **Structure**

These scripts require to have installed **pyseg/surf_dst** package, see:

> Martinez et al. (2018) Template-free detection and classification of heterogeneous membrane-bound complexes in cryo-electron tomograms. **bioRxiv**. doi: [https://doi.org/10.1101/413484](https://doi.org/10.1101/413484)

#### organization

This folder contains the scripts for computing the statistics that describe Apr2/3 complex distributions within different sets of tomograms (W - Wave and IT - Inner Territory). It also contains the input data with particles localization and shape but, due to space limitations, the tomograms with the segmentations used as masks for the volumes of interest are not included.
- **in** (folder): the input data with particles localization and shape.
- **ltomos_generator_{W,IT}_all.py**: scripts for preparing the input data
- **arp_uni_1st_{W,IT}_all.py**: scripts for computing the first order statistics
- **arp_uni_2nd_IT_all.py**: scripts for computing the second order statistics
- **gen_graphs.py**: script for generating the final graphs

#### angles

This folder contains the scripts to study the distribution of angles between two arbitrary vectors.
- **ang_dist_analytic.pdf**: theoretical analysis
- **angles.py**: set of utilities
- **run_ang_dist.py**: script for generating the experimental analysis
