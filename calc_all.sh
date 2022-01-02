#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate CLIPS

while getopts c:a:f:T:V: flag
do
    case "${flag}" in
        c) Ion1=${OPTARG};;
        a) Ion2=${OPTARG};;
        f) Solv=${OPTARG};;
        T) TEMPERATURE=${OPTARG};;
        V) NSOLV=${OPTARG};;
    esac
done

for i in `seq 1 10`
do
	bash calc_FE.sh -c $Ion1 -a $Ion2 -f $Solv -V $NSOLV -T $TEMPERATURE -i $i
	cp FE_"$Ion1"_"$Ion2"_"$Solv".pdf FE_$i.pdf
done

convert -background white -alpha remove -alpha off -verbose -loop 1 -delay 50 -dispose 3 -density 150 FE_{1..10}.pdf FE_"$Ion1"_"$Ion2"_"$Solv".gif
rm -rf FE_{1..10}.pdf \#*

