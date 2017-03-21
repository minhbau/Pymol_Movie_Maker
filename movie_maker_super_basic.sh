#!/bin/bash
#execute our lovely script with the correct pymol, python library support and passed commandline parameters
#echo $1 # path to pdb file (input)
#echo $2 # path pse file
#echo $3 # path pymol-script file

#super basic mode only uses 3 parameters, pdbfilename, path_pse_file and pymol script file

#used in our pymolscript as prefix for our script files, $galaxy is set in the startup script of the galaxy server
export MOVIEMAKERPATH="$galaxy""tools/customTools/movie_maker/"

#include library path to pymol in pythonpath, so python knows about the pymol 1.8.4 module
#include current directory in pythonpath, so scripts are available to import
export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${MOVIEMAKERPATH}:${PYTHONPATH}"
#used in pymolscript to save a seccond output
export POLAR_INTERACTION_FILENAME="polar_interaction_partners.txt"
export MOVIE_SCRIPT_FILENAME="movie_maker_basic_script.pml"


if [[ $# -eq 4 ]]
    then
        /home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker.py" --input "$1" > /home/webservices/philipp/movie_maker.log
    else
        (>&2 echo "'Super Basic mode' failed, wrong number of parameters, got "$#" expected 4")
fi


# move created pymol session from current directory to output directory
mv basic_movie.pse "$2"

# move textfile from pymol-script to destination directory
mv $POLAR_INTERACTION_FILENAME "$3"

# move pml script to destination directory, in combination with session allows easy modification of the movie
mv $MOVIE_SCRIPT_FILENAME "$4"