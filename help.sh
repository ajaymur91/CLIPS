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

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

