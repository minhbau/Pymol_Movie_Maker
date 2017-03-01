#!/bin/bash
#execute our lovely script with the correct pymol and passed commondaline parameters
#activate conda environment

#echo $galaxy
echo $1 # path to pdb file (input)                                                         
echo $2 # residue number
echo $3 # 3-letter-code
echo $4 # chain
echo $5

arr=$(echo $5 | sed -e 's/\(.*\)\/.*/\1\//')

echo ${arr}

#include library path to pymol in pythonpath, so python knows about the pymol 1.8.4 module
export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${PYTHONPATH}"
export MOVIE_MAKER_PATH= $galaxy"tools/customTools/movie_maker/"


#conda environment no available
#source activate special_pymol

#$1 is path to pdb file used for the visualization
#$2 is a required output path
#echo $1 #input file
#echo $2 #output file

#use the pymol 1.8.4 installation for our script
/home/webservices/philipp/special_pymol/pymol -c -u $MOVIE_MAKER_PATH"movie_maker_basic.py" -- $1 $2 $2 $3 $4 $5 $galaxy > /home/webservices/philipp/movie_maker_basic.log

#WHAT?
cd ${arr}

# move created pymol session from current directory to output directory
mv basic_movie_1.2.pse $5

#deactivate conda env
#source deactivate
