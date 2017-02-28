#!/bin/bash
#execute our lovely script with the correct pymol and passed commondaline parameters
#activate conda environment

echo $1 # path to ppdb file (input)                                                         
echo $2 # residue number
echo $3 # 3-letter-code
echo $4 # chain
echo $5

arr=$(echo $5 | sed -e 's/\(.*\)\/.*/\1\//')

echo ${arr}

source activate special_pymol

#$1 is path to pdb file used for the visualization
#$2 is a required output path
#echo $1 #input file
#echo $2 #output file

pymol -c -u $galaxy$"tools/customTools/movie_maker/movie_maker_basic.py" -- $1 $2 $2 $3 $4 $5 $galaxy > /home/webservices/philipp/movie_maker_basic.log

#WHAT?
cd ${arr}

# move created pymol session from current directory to output directory
mv basic_movie_1.2.pse $5

#deactivate conda env
source deactivate
