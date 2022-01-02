![](CLIPS.gif)

# External Requirements:
- git 
(Linux: https://git-scm.com/download/linux, MacOS: https://git-scm.com/download/mac)

- anaconda
(Linux: https://docs.anaconda.com/anaconda/install/linux/, MacOS: https://docs.anaconda.com/anaconda/install/mac-os/)

# Getting CLIPS code
- tested on linux-64 (Ubuntu 20.04.3, CentOS-7) and osx-64 (tested on macOS Catalina)
```
git clone ssh://git@github.com/ajaymur91/CLIPS.git
conda env create -f env.yaml # For linux-64
conda env create -f mac-env.yaml # For osx-64
```
# Usage (defaults)
```
cd CLIPS 
bash CLIPS.sh (uses defaults, results stored in FE.pdf, FE.gif, barrier, bindE)
```
# Usage (non-defaults)
  - Select some input parameters in CLIPS.sh
  - or use optional command line args:
```
	Syntax: bash CLIPS.sh -c LI -a TF -f EC -n 1 -T 313 -P LI -N S1 -S O2 -R 2 (results stored in FE_Ion1_Ion2_Solv.pdf, FE_Ion1_Ion2_Solv.gif, barrier, bindE)
	Options:												
     c) cation (LI or NA) (default = LI)									
     a) anion (TF or TFSI or BLB) (default = TF)								
     f) solvent (EC or SOL) (default is EC = Ethylene carbonate, SOL = water)					
     n) OpenMP processes (1 or 2 or 3 ...)									
     T) TEMPERATURE (default = 313)										
     P) cation reference atom for distance measurement (LI or NA) (default = LI)				
     N) anion reference atom (S1 for TF/TFSI, CL for Chloride) (default = S1)					
     S) solvent reference atom (O2 for EC, OW for water) (default = O2)					
     R) R_SOL (Confines solvent within sphere around center of ions (default = 2 nm)				
     h) print this help	

```
  - Currently supports LI, TFSI, TF, NA, CL, EC and water(SOL) only.
  - (Figshare animation) [https://doi.org/10.6084/m9.figshare.16755277.v3](https://doi.org/10.6084/m9.figshare.16755277.v3) (CLIPS.gif)
  - (Preprint) https://arxiv.org/abs/2110.14036 
  - (Citation) Muralidharan, Ajay, and Arun Yethiraj. "Fast estimation of ion-pairing for screening electrolytes: A cluster can approximate a bulk liquid." arXiv preprint arXiv:2110.14036 (2021).
