
#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate CLIPS

# Stop if you encounter error
set -e

Help()
{
   # Display Help
   echo "														"
   echo "	Syntax: bash NSOLV.sh -f EC -n 1 -T 313 								"
   echo "	Options:												"
   echo "     f) solvent (EC or SOL) (default is EC = Ethylene carbonate, SOL = water)					"
   echo "     n) OpenMP processes (1 or 2 or 3 ...)									"
 #  echo "     T) TEMPERATURE (default = 313)										"
   echo "     h) print this help											"		
   echo "														"
}

# Defaults
#TEMPERATURE=313
Solv=EC
NTOMP=1
Ns=1000       # Attempt to insert Ns solvents (for solvent density calculation) 
Time=0.05      # Equil. time in ns
NSTEPS=$(echo "$Time" | awk '{ printf("%.0f", $1*500000) }')
Rfill=0.90

# Inputs
while getopts f:n:h flag
do
    case "${flag}" in
        f) Solv=${OPTARG};;
        n) NTOMP=${OPTARG};;
        h) # display Help
           Help
           exit;;
    esac
done

# Print 
#echo "Solv: $Solv";
#echo "NTOMP: $NTOMP";
#echo "TEMPERATURE: $TEMPERATURE";
#echo "R_SOL: $R_SOL";
#echo "NSOLV: $NSOLV";

# Retry function  
function fail {
  echo $1 >&2
  exit 1
}

function retry {
  local n=1
  local max=3
  local delay=1
  while true; do
    "$@" && break || {
      if [[ $n -lt $max ]]; then
        ((n++))
	echo "CLIPS failed. Retry Attempt $((n+1))/$max:"
        sleep $delay;
      else
        fail "The command has failed after $max attempts."
      fi
    }
  done
}

# Insert Ns solvents
Ns=$(gmx insert-molecules -f struct/$Solv.gro -ci struct/$Solv.gro -o den.gro -box 2.1 2.1 2.1 -nmol 500 -try 200 -scale 0.57 2> /dev/stdout  | grep "Output configuration contains" | awk '{ print $(NF-1) }')

cat << EOF > density.top
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               3               yes             0.5     0.5

#include "itp/atomtypes.itp"
#include "itp/$Solv.itp"

[ system ]
; name
$Solv-density-calc

[ molecules ]
;name number
$Solv   $Ns
EOF

# Minimize and equilibrate neat solvent
gmx grompp -f mdp/min.mdp -c den.gro -p density.top -o density_min.tpr &> /dev/null
retry gmx mdrun -deffnm density_min -ntomp 1 -nsteps 50000 &> /dev/null
gmx grompp -f mdp/npt.mdp -c density_min.gro -p density.top -o density_eq.tpr &> /dev/null
retry gmx mdrun -deffnm density_eq -ntomp $NTOMP -nsteps $NSTEPS &> /dev/null

# Identify number of solvents required for CLIPS program ()
V=$(tail -n 1 density_eq.gro | awk '{print $1^3}')
NSOLV=$(echo $Ns*4*3.1416*0.3333*$Rfill^3/$V | bc)
echo $NSOLV | awk '{ printf("%.0f", $1) }'
