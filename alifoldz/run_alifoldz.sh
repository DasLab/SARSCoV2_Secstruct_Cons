#!/bin/bash

for i in $(seq 0 751); do
	echo $i
	alifoldz.pl --forward < "windows/alignment_"$i".fa" | tail -n 2 | head -n 1 | awk '{ print $4,$5,$6,$7 }' >> alifoldz_data.dat
done
