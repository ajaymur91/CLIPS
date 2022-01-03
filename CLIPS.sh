set -e
#!/bin/bash
Help()
{
   # Display Help
   echo "														"
   echo "	Syntax: bash CLIPS.sh -c LI -a TF -f EC -n 1 -T 313 -P LI -N S1 -S O2 -R 2 				"
   echo "	Options:												"
   echo "     c) cation (LI or NA) (default = LI)									"		
   echo "     a) anion (TF or TFSI or BLB) (default = TF)								"
   echo "     f) solvent (EC or SOL) (default is EC = Ethylene carbonate, SOL = water)					"
   echo "     n) OpenMP processes (1 or 2 or 3 ...)									"
   echo "     T) TEMPERATURE (default = 313)										"
   echo "     t) sampling time (default = 10 ns)										"
   echo "     P) cation reference atom for distance measurement (LI or NA) (default = LI)				"
   echo "     N) anion reference atom (S1 for TF/TFSI, CL for Chloride) (default = S1)					"
   echo "     S) solvent reference atom (O2 for EC, OW for water) (default = O2)					"
   echo "     R) R_SOL (Confines solvent within sphere around center of ions (default = 2 nm)				"
   echo "     h) print this help											"		
   echo "														"
}

# Defaults
TEMPERATURE=313
Ion1=LI
Ion2=TF
Solv=EC  
NTOMP=1   #OMP THREADS
R_SOL=2   #Solvent barrier
CA1=LI    #Solute_reference_Atom
CA2=S1    #Solute_reference_Atom
SA21=C4   #Solvent_reference_Atom 
time=10   #ns

while getopts c:a:f:n:T:t:P:N:S:R:h flag
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
        h) # display Help
           Help
           exit;;
    esac
done

# Print Inputs
echo "bash CLIPS.sh -h (Prints help function to set no defaults)";
echo "Cation: $Ion1";
echo "Anion: $Ion2";
echo "Solv: $Solv";
echo "NTOMP: $NTOMP";
echo "TEMPERATURE: $TEMPERATURE";
echo "R_SOL: $R_SOL";
echo " ";

# Estimate NSOLV
echo "Estimating number of solvent molecules required for CLIPS cluster calculation ";
NSOLV=$(bash NSOLV.sh -f $Solv -n $NTOMP)
echo "NSOLV: $NSOLV";
sleep 1

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
	echo "CLIPS failed. Retry Attempt $((n+1))/$max:"
        sleep $delay;
      else
        fail "The command has failed after $max attempts."
      fi
    }
  done
}

# Run CLIPS
time retry bash run_n.sh -c $Ion1 -a $Ion2 -f $Solv -n $NTOMP -T $TEMPERATURE -t $time -P $CA1 -N $CA2 -S $SA21 -R $R_SOL -V $NSOLV
