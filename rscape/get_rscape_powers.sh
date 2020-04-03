#!/bin/bash

cd windows/
for ii in $(seq 1 752); do
	echo $ii
	R-scape --fold "alignment_"$ii".sto"
	total_power=$(cat "alignment_"$ii"_1.fold.power" | grep -v "#" | tail -n +2 | awk '{ print $NF }' | paste -s -d+ | bc)
	num_bp=$(cat "alignment_"$ii"_1.fold.power" | grep -v "#" | tail -n +2 | awk '{ print $NF }' | wc -l)
	alignment_power=$(echo "scale=5; "$total_power"/"$num_bp | bc | awk '{printf "%.5f", $0}')
	echo -e $ii"\t"$total_power"\t"$num_bp"\t"$alignment_power >> alignment_powers.dat
	rm "alignment_"$ii"_1"*
done
