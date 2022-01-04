#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate CLIPS

# Stop if you encounter error
set -e
export GMX_MAXBACKUP=-1     # Overwrites
export PLUMED_MAXBACKUP=-1  # Unlimited backups
###############################

while getopts c:a:f:n:T:t:P:N:S:R:V: flag
do
    case "${flag}" in
        c) Ion1=${OPTARG};;
        a) Ion2=${OPTARG};;
        f) Solv=${OPTARG};;
        n) NTOMP=${OPTARG};;
        T) TEMPERATURE=${OPTARG};;
        t) time=${OPTARG};;
        P) CA1=${OPTARG};;
        N) CA2=${OPTARG};;
        S) SA21=${OPTARG};;
        R) R_SOL=${OPTARG};;
        V) NSOLV=${OPTARG};;
    esac
done

# Default Inputs 
#TEMPERATURE=313
AT=0.65
KAPPA=20000
#Ion1=LI
#Ion2=TF
#Solv=EC
#NTOMP=2
Ion_q1=1
Ion_q2=-1
Solv_q=0
#CA1=LI   #Solute_central_Atom
#CA2=S1   #Solute_central_Atom
#SA21=C4  #Solvent_binding_Atom 
#NSOLV=30
Tot_q=$(($Ion_q1+$Ion_q2))

#Inner shell radius (nm)
R1=0.28
R2=0.28
#R_SOL=2.0

#Trajectory sampling time (fs) - do not change
nsteps=50000        # Minimization
nstepsmd=100000     # Equilibration
nstepsmtd=$(echo "$time*500000" | bc | awk '{printf "%.0f", $1}')   # Enhanced Sampling steps

# Parse non-default inputs if available

echo "Cation: $Ion1";
echo "Anion: $Ion2";
echo "Solv: $Solv";
echo "NTOMP: $NTOMP";
echo "TEMPERATURE: $TEMPERATURE";
echo "sampling time: $time (ns)";
echo "R_SOL: $R_SOL";
echo "NSOLV: $NSOLV";
###############################

cat << EOF > ion.top
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5

#include "itp/atomtypes.itp"
#include "itp/$Ion1.itp"
#include "itp/$Ion2.itp"
#include "itp/$Solv.itp"

[ system ]
; name
$Ion1-$Ion2-$Solv-GP 

[ molecules ]
;name number
$Ion1      1
$Ion2      1
EOF

###############################

cat << EOF > verlet.mdp
define                  = -DFLEXIBLE
integrator              = sd        
nsteps                  = 1000000   
dt                      = 0.002     

nstxout                 = 2000     
nstvout                 = 2000     
nstenergy               = 2000     
nstcalcenergy           = 1
nstlog                  = 2000     

cutoff-scheme            = verlet
rlist                    = 1.5   
rcoulomb                 = 1.5   
rvdw                     = 1.5   
pbc                      = xyz
ns-type                  = simple
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = all-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2                   ; time constant, in ps
ref_t                   = 313                   ; reference temperature, one for each group, in K

EOF

###############################
cat << EOF > min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 0.0001          ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.00001          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
nstenergy                = 10
ns-type                  = simple
continuation             = no  ; does the same thing as unconstrained_start
EOF

###############################


# Start from Ion.gro and insert waters and generate force-field
cp ion.top system.top
gmx insert-molecules -f struct/"$Ion1".gro -ci struct/$Ion2.gro -o Ions.gro -box 1.5 1.5 1.5 -nmol 1 -try 1000 -scale 3 

NSOLV=$(gmx insert-molecules -f Ions.gro -ci struct/$Solv.gro -o IonW.gro -box 1.8 1.8 1.8 -nmol $NSOLV -try 1000 -scale 0.57 2> /dev/stdout  | grep "Output configuration contains" | awk '{ print $(NF-1)-2 }')
cat << EOF >> system.top
$Solv   $NSOLV
EOF

###############################

# Make index and plumed.dat (to create cluster)
gmx select -f IonW.gro -s IonW.gro -on CA1.ndx -select "atomname $CA1 and resnr 1" 
gmx select -f IonW.gro -s IonW.gro -on CA2.ndx -select "atomname $CA2 and resnr 2" 
gmx select -f IonW.gro -s IonW.gro -on mtd.ndx -select "atomname $SA21 and resnr > 2" 

Nt=$NSOLV 
cat << EOF > plumed.dat
CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
com: CENTER ATOMS=CA1,CA2
SA: GROUP NDX_FILE=mtd.ndx # NDX_GROUP=atomname__and_resnr_<_6
cn: COORDINATION GROUPA=com GROUPB=SA R_0=$R_SOL NN=12
LW: LOWER_WALLS ARG=cn AT=${Nt} KAPPA=1000

di: DISTANCE ATOMS=CA1,CA2
UPPER_WALLS ...
 ARG=di
 AT=${AT}
 KAPPA=${KAPPA}
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS
RESTRAINT ARG=di AT=0.6 KAPPA=5000.0 LABEL=restraint
EOF
echo "PRINT ARG=* FILE=COLVAR STRIDE=100" >> plumed.dat

###############################

# Make box bigger
	gmx editconf -f IonW.gro -c -box 4 4 4 -o start1.gro 
# Minimize
	echo -e "\n Run Minimization $Ion1-$Ion2 and $NSOLV $Solv \n" 
	gmx grompp -f min.mdp -c start1.gro -p system.top -o min.tpr 
	gmx mdrun -deffnm min -nsteps 100000 -plumed plumed.dat 
# Make box bigger
	gmx editconf -f min.gro -c -box 4.97 4.97 4.97 -o start.gro 
# Equilibrate
	echo -e "\n Run MD - $Ion1-$Ion2 - $NSOLV $Solv \n" 
	gmx grompp -f verlet.mdp -c start.gro -p system.top -o md.tpr 
	gmx mdrun -deffnm md -nsteps $nstepsmd -plumed plumed.dat -ntomp $NTOMP

###############################

# Create plumed input for running OPES
cat << EOF > plumed_MTD.dat
CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
com: CENTER ATOMS=CA1,CA2
SA: GROUP NDX_FILE=mtd.ndx # NDX_GROUP=atomname__and_resnr_<_6
cn: COORDINATION GROUPA=com GROUPB=SA R_0=$R_SOL NN=12
EOF

echo "LOWER_WALLS ARG=cn AT=${Nt} KAPPA=10 LABEL=LW" >> plumed_MTD.dat

cat << EOF >> plumed_MTD.dat
di: DISTANCE ATOMS=CA1,CA2

opes: OPES_METAD ...
  ARG=di
  FILE=Kernels.data
  TEMP=${TEMPERATURE}
  PACE=50
  BARRIER=100
  #SIGMA=0.05
  #SIGMA_MIN=0.0005
  STATE_WFILE=State.data
  STATE_WSTRIDE=10000
  STORE_STATES
 ...

UPPER_WALLS ...
 ARG=di
 AT=${AT}
 KAPPA=${KAPPA}
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS

PRINT FMT=%g STRIDE=10 FILE=Colvar.data ARG=di,cn,*.bias
EOF

###############################

# WT-MTD (1000 ps)
echo -e "\n Run WT-MTD - $Ion1 - $Ion2 and $NSOLV $Solv \n"
gmx grompp -f verlet.mdp -c md.gro -p system.top -o mtd.tpr 
gmx mdrun -deffnm mtd -nsteps $nstepsmtd -plumed plumed_MTD.dat -ntomp $NTOMP 

###############################
rm -rf barrier \#*
bash calc_all.sh -c $Ion1 -a $Ion2 -f $Solv -T $TEMPERATURE -V $NSOLV

