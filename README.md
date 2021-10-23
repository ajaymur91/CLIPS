# Calculate dissociation barrier using OPES

  ```
# Requirements:
conda env create -f env.yaml
  ```
  ```
# Usage
conda activate CLIPS
bash run_n.sh
bash run_n.sh -c LI -a TFSI -f EC -n 2 -T 313 -P LI -N S1 -S C4
  ```
  - Select some input parameters in run_n.sh or
  -  optional command line args:
  -     c) cation (LI or NA)
        a) anion (TF or TFSI or BLB)
        f) solvent (EC or SOL) (SOL stands for water)
        n) OpenMP processes (1 or 2 or 3 ...) 
        T) TEMPERATURE (313 default)
        P) (cation reference atom for distance measurement)
        N) (anion reference atom for distance measurement)
        S) (solvent reference atom for distance measurement)

  - Currently supports LI, TFSI, OTF, NA, CL, EC and water only.
  - (Figshare animation) https://doi.org/10.6084/m9.figshare.16755277.v3 (CLIPS.gif)

   ![](CLIPS.gif)

  ```
