   ![](CLIPS.gif)
# Calculate dissociation barrier using OPES

  ```
# External Requirements:
- git 
- anaconda
- tested on linux-64 and osx-64
  ```
  ```
# Getting CLIPS code
git clone ssh://git@github.com/ajaymur91/CLIPS.git
conda env create -f env.yaml # For linux-64
conda env create -f mac-env.yaml # For osx-64

# Usage
cd CLIPS 
bash run_n.sh (uses defaults, results stored in FE.pdf)
  - Select some input parameters in run_n.sh
  - or use optional command line args:

bash run_n.sh -c LI -a TFSI -f EC -n 2 -T 313 -P LI -N S1 -S C4 (setting non defaults)

        c) cation (LI or NA)
        a) anion (TF or TFSI or BLB)
        f) solvent (EC or SOL) (SOL stands for water)
        n) OpenMP processes (1 or 2 or 3 ...) 
        T) TEMPERATURE (313 default)
        P) cation reference atom for distance measurement (LI or NA)
        N) anion reference atom for distance measurement (S1 for TF/TFSI, CL for Chloride)
        S) solvent reference atom for distance measurement (O2 for EC, OW for water)
        R) R_SOL (Confines solvent within sphere around ions, default=1.8nm)
        V) N_SOLV (currently empirical: 30 for EC and 80 for water. Anything greater than 20% of bulk density works.)

  - Currently supports LI, TFSI, OTF, NA, CL, EC and water only.
  - (Figshare animation) https://doi.org/10.6084/m9.figshare.16755277.v3 (CLIPS.gif)

  ```
