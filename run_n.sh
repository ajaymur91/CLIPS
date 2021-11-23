# Script for dimer interactions
# Requires anaconda
# See ./readme_conda 

#!/bin/bash
# Stop if you encounter error
eval "$(conda shell.bash hook)"
conda activate CLIPS
set -e
export GMX_MAXBACKUP=-1     # Overwrites
export PLUMED_MAXBACKUP=-1  # Unlimited backups
#conda activate nucl
#export OMP_NUM_THREADS=6
###############################

#source gmx_plumed/bin/activate 
# Edit inputs (defaults) below
TEMPERATURE=313
Ion1=LI
Ion2=TF
Solv=EC
NTOMP=3
Ion_q1=1
Ion_q2=-1
Solv_q=0
CA1=LI   #Solute_central_Atom
CA2=S1   #Solute_central_Atom
SA21=C4  #Solvent_binding_Atom 
NSOLV=30
Tot_q=$(($Ion_q1+$Ion_q2))

#Inner shell radius (nm)
R1=0.28
R2=0.28
R_SOL=2

#Trajectory sampling time (ps) - do not change
   nsteps=50000
 nstepsmd=500000     #dt is set to 0.5 fs in mdp file
nstepsmtd=5000000 #0

while getopts c:a:f:n:T:P:N:S: flag
do
    case "${flag}" in
        c) Ion1=${OPTARG};;
        a) Ion2=${OPTARG};;
        f) Solv=${OPTARG};;
        n) NTOMP=${OPTARG};;
        T) TEMPERATURE=${OPTARG};;
        P) CA1=${OPTARG};;
        N) CA2=${OPTARG};;
        S) SA21=${OPTARG};;
    esac
done
echo "Cation: $Ion1";
echo "Anion: $Ion2";
echo "Solv: $Solv";
echo "NTOMP: $NTOMP";
echo "TEMPERATURE: $TEMPERATURE";

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

cat << EOF > md.mdp
define                  = -DFLEXIBLE
integrator              = sd        ; leap-frog integrator
nsteps                  = 1000000     ; 2 * 50000 = 100 ps
dt                      = 0.001     ; 0.5 fs

; Output control
nstxout                 = 2000     ; save coordinates every 1.0 ps
nstvout                 = 2000     ; save velocities every 1.0 ps
nstenergy               = 2000     ; save energies every 1.0 ps
nstcalcenergy           = 1
nstlog                  = 2000       ; update log file every 1.0 ps


; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
nstlist                  = 0
ns-type                  = simple
constraints              = none      ; bonds involving H are constrained
continuation             = no  ; does the same thing as unconstrained_start

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).

; Run parameters

; Temperature coupling is on
;tcoupl                  = nose-hoover           ; modified Berendsen thermostat
tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2                   ; time constant, in ps
ref_t                   = 313                   ; reference temperature, one for each group, in K

EOF
###############################


cat << EOF > verlet.mdp
define                  = -DFLEXIBLE
integrator              = sd        ; leap-frog integrator
nsteps                  = 1000000     ; 2 * 50000 = 100 ps
dt                      = 0.001     ; 0.5 fs

; Output control
nstxout                 = 2000     ; save coordinates every 1.0 ps
nstvout                 = 2000     ; save velocities every 1.0 ps
nstenergy               = 2000     ; save energies every 1.0 ps
nstcalcenergy           = 1
nstlog                  = 2000       ; update log file every 1.0 ps


; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = verlet    ; Buffered neighbor searching
rlist                    = 1.2   ;gas ph min (0 means no cutoff)
rcoulomb                 = 1.2   ;gas ph min (0 means no cutoff)
rvdw                     = 1.2   ;gas ph min (0 means no cutoff)
pbc                      = xyz
;nstlist                  = 10
ns-type                  = simple
;constraints              = none      ; bonds involving H are constrained
;continuation             = no  ; does the same thing as unconstrained_start
;DispCorr                = EnerPres
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = all-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).

; Run parameters

; Temperature coupling is on
;tcoupl                  = nose-hoover           ; modified Berendsen thermostat
tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2                   ; time constant, in ps
ref_t                   = 313                   ; reference temperature, one for each group, in K

EOF
###############################

cat << EOF > min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 0.1          ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
;cutoff-scheme            = verlet    ; Buffered neighbor searching
;rlist                    = 1.2   ;gas ph min (0 means no cutoff)
;rcoulomb                 = 1.2   ;gas ph min (0 means no cutoff)
;rvdw                     = 1.2   ;gas ph min (0 means no cutoff)
;pbc                      = xyz
;nstlist                  = 10
;ns-type                  = simple
;constraints              = none      ; bonds involving H are constrained
;continuation             = no  ; does the same thing as unconstrained_start

cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
;nstenergy                = 10
;nstxout                  = 10
;nstlist                  = 10
ns-type                  = simple
continuation             = no  ; does the same thing as unconstrained_start
;unconstrained_start     = yes ; depricated

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).
EOF

###############################

# Start from Ion.gro and insert waters and generate force-field
cp ion.top system.top
gmx insert-molecules -f struct/"$Ion1".gro -ci struct/$Ion2.gro -o Ions.gro -box 1.5 1.5 1.5 -nmol 1 -try 1000 -scale 3 #&> /dev/null
#cat << EOF >> system.top
#$Ion2   1
#EOF

gmx insert-molecules -f Ions.gro -ci struct/$Solv.gro -o IonW.gro -box 2 2 2 -nmol $NSOLV -try 1000 -scale 0.1 #&> /dev/null
cat << EOF >> system.top
$Solv   $NSOLV
EOF

###############################

# Make index and plumed.dat ( to enforce QCT criterion )
gmx select -f IonW.gro -s IonW.gro -on CA1.ndx -select "atomname $CA1 and resnr 1" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on CA2.ndx -select "atomname $CA2 and resnr 2" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on mtd.ndx -select "atomname $SA21 and resnr > 2" &> /dev/null

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
 AT=0.35
 KAPPA=2000.0
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS
EOF
echo "PRINT ARG=* FILE=COLVAR STRIDE=100" >> plumed.dat

###############################

# Minimize
gmx editconf -f IonW.gro -c -box 4 4 4 -o start1.gro #&> /dev/null
echo -e "\n Run Minimization $Ion1-$Ion2 and $NSOLV $Solv \n"
gmx grompp -f min.mdp -c start1.gro -p system.top -o min.tpr  #&> /dev/null
gmx mdrun -v -deffnm min -nsteps 10000 -plumed plumed.dat  #&> /dev/null

# Make box bigger (does not really matter - but do it anyway to be safe)
gmx editconf -f min.gro -c -box 4.786 4.786 4.786 -o start.gro #&> /dev/null

#Md (100 ps)
echo -e "\n Run MD - $Ion1-$Ion2 - $NSOLV $Solv \n"
 gmx grompp -f verlet.mdp -c start.gro -p system.top -o md.tpr #&> /dev/null
 gmx mdrun -v -deffnm md -nsteps $nstepsmd -plumed plumed.dat -ntomp $NTOMP #&> /dev/null

###############################
Nt=$NSOLV 

cat << EOF > plumed_MTD.dat
ene: ENERGY
CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
com: CENTER ATOMS=CA1,CA2
SA: GROUP NDX_FILE=mtd.ndx # NDX_GROUP=atomname__and_resnr_<_6
cn: COORDINATION GROUPA=com GROUPB=SA R_0=$R_SOL NN=12
EOF

echo "LOWER_WALLS ARG=cn AT=${Nt} KAPPA=1000 LABEL=LW" >> plumed_MTD.dat

cat << EOF >> plumed_MTD.dat
di: DISTANCE ATOMS=CA1,CA2

opes: OPES_METAD ...
  ARG=di
  FILE=Kernels.data
  TEMP=${TEMPERATURE}
  PACE=50
  BARRIER=100
  SIGMA=0.05
  SIGMA_MIN=0.0001
  #BIASFACTOR=25
  STATE_WFILE=State.data
  STATE_WSTRIDE=10000
  STORE_STATES
 ...

UPPER_WALLS ...
 ARG=di
 AT=0.65
 KAPPA=200000.0
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS

PRINT FMT=%g STRIDE=10 FILE=Colvar.data ARG=di,ene,cn,*.bias,opes.*
EOF

###############################

# WT-MTD (1000 ps)
echo -e "\n Run WT-MTD - $Ion1 - $Ion2 and $NSOLV $Solv \n"
 gmx grompp -f verlet.mdp -c md.gro -p system.top -o mtd.tpr -t md.cpt #&> /dev/null
 gmx mdrun -v -deffnm mtd -nsteps $nstepsmtd -plumed plumed_MTD.dat -ntomp $NTOMP #&> /dev/null

###############################
rm -rf barrier
bash calc_FE.sh $Ion1 $Ion2; 

# optional (view results)
#conda activate MUPDF
#mupdf-gl FE.pdf
