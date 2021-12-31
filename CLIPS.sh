
#!/bin/bash
Help()
{
   # Display Help
   echo "														"
   echo "	Syntax: bash run_n.sh -c LI -a TF -f EC -n 1 -T 313 -P LI -N S1 -S O2 -R 2 -V 30			"
   echo "	Options:												"
   echo "     c) cation (LI or NA) (default = LI)									"		
   echo "     a) anion (TF or TFSI or BLB) (default = TF)								"
   echo "     f) solvent (EC or SOL) (default is EC = Ethylene carbonate, SOL = water)					"
   echo "     n) OpenMP processes (1 or 2 or 3 ...)									"
   echo "     T) TEMPERATURE (default = 313)										"
   echo "     P) cation reference atom for distance measurement (LI or NA) (default = LI)				"
   echo "     N) anion reference atom (S1 for TF/TFSI, CL for Chloride) (default = S1)					"
   echo "     S) solvent reference atom (O2 for EC, OW for water) (default = O2)					"
   echo "     R) R_SOL (Confines solvent within sphere around center of ions (default = 2 nm)				"
   echo "     V) N_SOLV (currently empirical: 30 for EC and 80 for water. (greater than 20% of bulk density works.)	"		
   echo "     h) print this help											"		
   echo "														"
}

# Defaults
TEMPERATURE=313
AT=0.65
KAPPA=1000
Ion1=LI
Ion2=TF
Solv=EC
NTOMP=2
CA1=LI    #Solute_reference_Atom
CA2=S1    #Solute_reference_Atom
SA21=C4   #Solvent_reference_Atom 
NSOLV=30

while getopts c:a:f:n:T:P:N:S:R:V:h flag
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
        R) R_SOL=${OPTARG};;
        V) NSOLV=${OPTARG};;
        h) # display Help
           Help
           exit;;
    esac
done

# Print Inputs
echo "Cation: $Ion1";
echo "Anion: $Ion2";
echo "Solv: $Solv";
echo "NTOMP: $NTOMP";
echo "TEMPERATURE: $TEMPERATURE";
echo "R_SOL: $R_SOL";
echo "NSOLV: $NSOLV";

# Retry function  
function fail {
  echo $1 >&2
  exit 1
}

function retry {
  local n=1
  local max=10
  local delay=1
  while true; do
    "$@" && break || {
      if [[ $n -lt $max ]]; then
        ((n++))
        echo "Command failed. Attempt $n/$max:"
        sleep $delay;
      else
        fail "The command has failed after $n attempts."
      fi
    }
  done
}


#time retry bash run_n.sh -c $Ion1 -a $Ion2 -f $Solv -n $NTOMP -P $CA1 -N $CA2 -S $SA1 -R $R_SOL -V $NSOLV
