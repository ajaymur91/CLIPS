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
  -     c) Ion1=${OPTARG};; cation (LI or NA)
        a) Ion2=${OPTARG};; anion (TF or TFSI or BLB)
        f) Solv=${OPTARG};; solvent (EC or SOL) (SOL stands for water)
        n) NTOMP=${OPTARG};; OpenMP processes (1 or 2 or 3 ...) 
        T) TEMPERATURE=${OPTARG};; (313 default)
        P) CA1=${OPTARG};; (cation reference atom for distance measurement)
        N) CA2=${OPTARG};; (anion reference atom for distance measurement)
        S) SA21=${OPTARG};; (solvent reference atom for distance measurement)

  - Currently supports LI, TFSI, OTF, NA, CL, EC and water only.
  - (Figshare animation) https://doi.org/10.6084/m9.figshare.16755277.v3 (CLIPS.gif)

   ![](CLIPS.gif)

  ```
