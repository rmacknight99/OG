#!/bin/bash
# Script to get matching complexes 


for file in data/experimental_spectras/*.csv; do
	fullname="$(basename -- $file)"
	basename=${fullname%.csv}
	python src/get_matching.py $basename

done


