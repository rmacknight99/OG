#!/bin/bash


mkdir ../experimental_spectras
find ../Discovery_EDA/EDA/UV/ -type f -name *.Sample.Raw.asc -exec cp {} ../experimental_spectras \;
	

