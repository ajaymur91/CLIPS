# Script for dimer interactions
# Requires anaconda
# See ./readme_conda 

#!/bin/bash
# Stop if you encounter error
set -e
export GMX_MAXBACKUP=-1     # Overwrites
export PLUMED_MAXBACKUP=-1  # Unlimited backups
#export OMP_NUM_THREADS=6
###############################

#source gmx_plumed/bin/activate 
# Edit inputs below
TEMPERATURE=313
POFF=${1:-0}
NTOMP=${2:-3}

#Ion1 
Ion1=LI
Ion_q1=1
Ion2=TFSI
Ion_q2=-1

#Solvent
Solv=EC
Solv_q=0

Tot_q=$(($Ion_q1+$Ion_q2))

#Solute_central_Atom
CA1=LI
CA2=S1

#Solvent_binding_Atom 
SA11=O2
SA12=O2
SA21=C4
SA22=C4
SA23=C4
SA24=C4

#Number of Solv
n1=20
n2=20

#Inner shell radius (nm)
R1=0.28
R2=0.28
R_SOL=2
#FUNCTIONAL and BASIS SET
#FUNCTIONAL=MP2
#BASIS_SET=aug-cc-pvdz

#Trajectory sampling time (ps) - do not change
nsteps=50000
nstepsmd=100000     #dt is set to 0.5 fs in mdp file
nstepsmtd=10000000 #0
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
rlist                    = 1.5   ;gas ph min (0 means no cutoff)
rcoulomb                 = 1.5   ;gas ph min (0 means no cutoff)
rvdw                     = 1.5   ;gas ph min (0 means no cutoff)
pbc                      = xyz
nstlist                  = 10
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

cat << EOF > min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 0.001          ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.001          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
;cutoff-scheme            = verlet    ; Buffered neighbor searching
;rlist                    = 1.5   ;gas ph min (0 means no cutoff)
;rcoulomb                 = 1.5   ;gas ph min (0 means no cutoff)
;rvdw                     = 1.5   ;gas ph min (0 means no cutoff)
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
nstenergy                = 10
nstxout                  = 10
nstlist                  = 0
ns-type                  = simple
continuation             = no  ; does the same thing as unconstrained_start
;unconstrained_start     = yes ; depricated

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).
EOF

###############################

# Start from Ion.gro and insert waters and generate force-field
cp ion.top system.top
gmx insert-molecules -f struct/"$Ion1".gro -ci struct/$Ion2.gro -o Ions.gro -box 1 1 1 -nmol 1 -try 1000 -scale 3 #&> /dev/null
#cat << EOF >> system.top
#$Ion2   1
#EOF

gmx insert-molecules -f Ions.gro -ci struct/$Solv.gro -o IonW.gro -box 2 2 2 -nmol $((n1+n2)) -try 1000 -scale 0.1 #&> /dev/null
cat << EOF >> system.top
$Solv   $((n1+n2))
EOF

###############################

# Make index and plumed.dat ( to enforce QCT criterion )
gmx select -f IonW.gro -s IonW.gro -on CA1.ndx -select "atomname $CA1 and resnr 1" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on CA2.ndx -select "atomname $CA2 and resnr 2" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on SA1.ndx -select "(atomname $SA11 || atomname $SA12) and resnr < $((3+$n1))" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on SA2.ndx -select "(atomname $SA21 || atomname $SA22 || atomname $SA23 || atomname $SA24 ) and resnr > $((2+$n1))" &> /dev/null
gmx select -f IonW.gro -s IonW.gro -on mtd.ndx -select "atomname $SA21 and resnr > 2" &> /dev/null

#cat << EOF > plumed.dat
#CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
#SA1: GROUP NDX_FILE=SA1.ndx # NDX_GROUP=atomname_${SA11}_||_and_resnr_<_$((3+$n1))
#cn1: COORDINATION GROUPA=CA1 GROUPB=SA1 R_0=$(echo $R1+0.12 | bc) NN=12
#CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
#SA2: GROUP NDX_FILE=SA2.ndx # NDX_GROUP=atomname_${SA2}_and_resnr_>_$((2+$n1))
#cn2: COORDINATION GROUPA=CA2 GROUPB=SA2 R_0=$(echo $R2+0.12 | bc) NN=12
#EOF
#
#echo "UPPER_WALLS ARG=cn1 AT=6.1 KAPPA=500 LABEL=UWcn1" >> plumed.dat
#echo "UPPER_WALLS ARG=cn2 AT=8.1 KAPPA=500 LABEL=UWcn2" >> plumed.dat
#echo "LOWER_WALLS ARG=cn1 AT=4.9 KAPPA=500 LABEL=LWcn1" >> plumed.dat
#echo "LOWER_WALLS ARG=cn2 AT=5.9 KAPPA=500 LABEL=LWcn2" >> plumed.dat

cat << EOF > plumed.dat
CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
com: CENTER ATOMS=CA1,CA2
SA: GROUP NDX_FILE=mtd.ndx # NDX_GROUP=atomname__and_resnr_<_6
cn: COORDINATION GROUPA=com GROUPB=SA R_0=1 NN=12

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
#k=0;
#for i in $(tail -n 1 SA1.ndx)
#do
#k=$(($k+1))
#echo "d$k: DISTANCE ATOMS=CA1,$i" >> plumed.dat
#echo -e "UPPER_WALLS ARG=d$k AT=$R1 KAPPA=500 LABEL=w$k" >> plumed.dat
#cat << EOF >> plumed.dat
#EOF
#done
#
#for i in $(tail -n 1 SA2.ndx)
#do
#k=$(($k+1))
#echo "d$k: DISTANCE ATOMS=CA2,$i" >> plumed.dat
#echo -e "UPPER_WALLS ARG=d$k AT=$R2 KAPPA=500 LABEL=w$k" >> plumed.dat
#cat << EOF >> plumed.dat
#EOF
#done
echo "PRINT ARG=* FILE=COLVAR STRIDE=100" >> plumed.dat
#echo "PRINT ARG=*.bias FILE=WALL STRIDE=100" >> plumed.dat

###############################

# Minimize
echo -e "\n Run Minimization $Ion1-$Ion2 and $((n1+n2)) $Solv \n"
gmx grompp -f min.mdp -c IonW.gro -p system.top -o min.tpr  #&> /dev/null
gmx mdrun -deffnm min -nsteps 1000000 -plumed plumed.dat  #&> /dev/null

# Make box bigger (does not really matter - but do it anyway to be safe)
gmx editconf -f min.gro -c -box 4 4 4 -o start.gro #&> /dev/null

#Md (100 ps)
echo -e "\n Run MD - $Ion1-$Ion2 - $((n1+n2)) $Solv \n"
 gmx grompp -f verlet.mdp -c start.gro -p system.top -o md.tpr #&> /dev/null
 gmx mdrun -v -deffnm md -nsteps $nstepsmd -plumed plumed.dat -ntomp $NTOMP #&> /dev/null

#gmx solvate -cp min2.gro -cs struct/ec_equil_1000.gro -p system.top -o solv.gro 
#gmx select -f solv.gro -s solv.gro -on solv.ndx -select "atomname O2" &> /dev/null
###############################
Nt=$((n1+n2)) 

cat << EOF > plumed_MTD.dat
CA1: GROUP NDX_FILE=CA1.ndx # NDX_GROUP=atomname_${CA1}_and_resnr_1
CA2: GROUP NDX_FILE=CA2.ndx # NDX_GROUP=atomname_${CA2}_and_resnr_2
com: CENTER ATOMS=CA1,CA2
SA: GROUP NDX_FILE=mtd.ndx # NDX_GROUP=atomname__and_resnr_<_6
cn: COORDINATION GROUPA=com GROUPB=SA R_0=$R_SOL NN=12
EOF

echo "LOWER_WALLS ARG=cn AT=${Nt} KAPPA=500 LABEL=LW" >> plumed_MTD.dat
#k=0;
#for i in $(tail -n 1 SA1.ndx)
#do
#k=$(($k+1))
#echo "d$k: DISTANCE ATOMS=CA1,$i" >> plumed_MTD.dat
#echo -e "UPPER_WALLS ARG=d$k AT=$R1 KAPPA=500 LABEL=w$k" >> plumed_MTD.dat
#done
#
#for i in $(tail -n 1 SA2.ndx)
#do
#k=$(($k+1))
#echo "d$k: DISTANCE ATOMS=CA2,$i" >> plumed_MTD.dat
#echo -e "UPPER_WALLS ARG=d$k AT=$R2 KAPPA=500 LABEL=w$k" >> plumed_MTD.dat
#done 

cat << EOF >> plumed_MTD.dat
di: DISTANCE ATOMS=CA1,CA2

#METAD ...
# LABEL=metad
# ARG=di
# SIGMA=0.08
# HEIGHT=1.
# BIASFACTOR=15
# TEMP=${TEMPERATURE}
# PACE=500
# GRID_MIN=0.0
# GRID_MAX=1.2
# GRID_BIN=500
# CALC_RCT
# #RCT_USTRIDE=50
# #REWEIGHTING_NGRID=500
#... METAD

opes: OPES_METAD ...
  ARG=di
  FILE=Kernels.data
  TEMP=${TEMPERATURE}
  PACE=100
  BARRIER=20
  SIGMA=0.04
  SIGMA_MIN=0.02
  #BIASFACTOR=25
  STATE_WFILE=State.data
  STATE_WSTRIDE=10000
  STORE_STATES
...



UPPER_WALLS ...
 ARG=di
 AT=0.65
 KAPPA=2000.0
 EXP=2
 EPS=1
 OFFSET=0.
 LABEL=uwall
... UPPER_WALLS

PRINT FMT=%g STRIDE=100 FILE=Colvar.data ARG=di,cn,*.bias,opes.*
# PRINT ARG=di,cn,*.bias,metad.rbias FILE=BIAS STRIDE=100
EOF

###############################

# WT-MTD (1000 ps)
echo -e "\n Run WT-MTD - $Ion1 - $Ion2 and $((n1+n2)) $Solv \n"
 gmx grompp -f verlet.mdp -c md.gro -p system.top -o mtd.tpr #&> /dev/null
 gmx mdrun -v -deffnm mtd -nsteps $nstepsmtd -plumed plumed_MTD.dat -ntomp $NTOMP #&> /dev/null

###############################
rm -rf barrier
bash calc_FE.sh; xdg-open FE.pdf
 
