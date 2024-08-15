#!/bin/bash

# This file contains a set of instructions to run the main file OPO.cpp.
# Use this file to perform simulations in bulk. This means, for example, 
# systematically varying the input power, the cavity reflectivity, etc.

TIMEFORMAT='It took %0R seconds.' 
time { 

clear      # Clear screen
rm *.dat   # This removes all .dat files in the current folder. Comment this line for safe.
rm OPO   # This removes a previuos executable file (if it exist)


#####################################################################
# -----Compilation-----
g++ OPO.cpp -o OPO -lfftw3 -lfftw3f -lm
FOLDERSIM="Simulations"


# -lfftw3 -lfftw3f        : flags needed to perform FFT
#####################################################################

# The variables defined below (ARGX) will be passed as arguments to the main file 
# OPO.cpp on each execution via the argv[X] instruction.

ARG1=(30)     	# Pump power                (ARG1)
ARG2=(29)      	# Beam waist (um)           (ARG2)
ARG3=(55.66) 	# 1st peltier temperature   (ARG3)

SOLS=1			# Save Only Last Slice (FALSE=0,TRUE=1)
STMP=0			# Save the entire temperature tensor (FALSE=0,TRUE=1)

x=$(awk 'BEGIN{for(i=20; i<=40; i=i+1)print i}')
# n=$(awk 'BEGIN{for(i=1; i<=31; i=i+3)print i}')
# t1=$(awk 'BEGIN{for(i=54; i<=59; i=i+0.1)print i}')
# Each for-loop span over one or more values defined in the previous arguments. 
# 		for N in $x
# 		do


for (( p=0; p<${#ARG1[@]}; p++ ))
do
# for N in $n
# do
	for (( x=0; x<${#ARG2[@]}; x++ ))
	do
# 	for WA in $x
# 	do
		for (( t1=0; t1<${#ARG3[@]}; t1++ ))
		do
# 		for TOVEN1 in $t1
# 		do

			N=${ARG1[$p]}
			printf "\nPower			= ${N} W\n" 
			
			WA=${ARG2[$x]}
			printf "\nBeam waist	= ${WA} \n" 
			
			TOVEN1=${ARG3[$t1]}
			printf "\nT peltier 1		= ${TOVEN1} ºC\n"
			
			printf "\nMaking directory...\n"
			FOLDER="sPPLT_P_${N}_W_T1_${TOVEN1}_WAIST_${WA}"
			FILE="sPPLT_P_${N}_W_T1_${TOVEN1}_WAIST_${WA}.txt"
			
			printf "Bash execution and writing output file...\n\n"
			./OPO $N $WA $TOVEN1 $SOLS $STMP| tee -a $FILE
			
			printf "Bash finished!!\n\n" 
			mkdir $FOLDER
			mv *.dat $FOLDER"/"
			mv FILE.txt $FOLDER"/"

		done
	done
done


if [ -d "$FOLDERSIM" ]; then
	echo "Moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
else

	mkdir $FOLDERSIM
	echo "Creating and moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
fi

}
