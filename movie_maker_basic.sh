#!/bin/bash
#execute our lovely script with the correct pymol and passed commondaline parameters
#echo $1 # path to pdb file (input)
#echo $2 # ligand name
#echo $3 # chain
#echo $4 # use colorblind friendly coloring?
#echo $5 # path pse file
#echo $6 # path pymol-script file



#used in our pymolscript as prefix for our script files
export MOVIEMAKERPATH="$galaxy""tools/customTools/movie_maker/"

#include library path to pymol in pythonpath, so python knows about the pymol 1.8.4 module
#include current directory in pythonpath, so scripts are available to import
export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${MOVIEMAKERPATH}:${PYTHONPATH}"
#used in pymolscript to save a seccond output
export POLAR_INTERACTION_FILENAME="polar_interaction_partners.txt"

echo $POLAR_INTERACTION_FILENAME

arr=$(echo $5 | sed -e 's/\(.*\)\/.*/\1\//')

#path+="/xbscore_output_pdf.py"
echo ${arr}

#echo "working on $MOVIEMAKERPATH"
#echo "executing $MOVIEMAKERPATH""movie_maker_basic.py"
/home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker_basic.py" -- $1 $2 $2 $3 $4 $5 $galaxy > /home/webservices/philipp/movie_maker_basic.log

# move created pymol session from current directory to output directory
mv basic_movie.pse $5
# move textfile from pymol-script to destination directory
mv $POLAR_INTERACTION_FILENAME $6
