#!/bin/bash



cd data/finger_prints
for dimer in *; do
	ID=$(basename "${dimer%.*}")
	#xyz=$ID
	#xyz+='.xyz'
	mol=$ID
	mol+='.mol'
	#echo $xyz
	#obabel $xyz -O $mol
	#mv $mol ../finger_prints
	python ../../src/finger_print.py --mol $mol
done
