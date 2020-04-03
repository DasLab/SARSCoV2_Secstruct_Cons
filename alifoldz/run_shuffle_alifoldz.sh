#!/bin/bash

for i in $(seq 0 751); do
	echo $i
	cat "windows/alignment_"$i".fa" | shuffle-aln.pl | alifoldz.pl --forward | tail -n 2 | head -n 1 | awk '{ print $4,$5,$6,$7 }' >> alifoldz_shuffle_data.dat
done
