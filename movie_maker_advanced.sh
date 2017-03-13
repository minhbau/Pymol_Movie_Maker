#!/bin/bash
#execute our lovely script with the correct pymol, python library support and passed commandline parameters
#echo $1  # path to pdb file (input)
#echo $2  # ligand name
#echo $3  # chain
#echo $4  # use colorblind friendly coloring?
#echo $5  # path pse file
#echo $6  # path pymol-script file
#echo $7  # binding_site_radius
#echo $8  # check_halogen_interaction
#echo $9  # water_in_binding_site
#echo $10 # color_carbon
#echo $11 # session_export_version
#
#echo $7
#echo $8
#echo $9
#echo $10
#echo $11
#echo $12
#echo $13

#used in our pymolscript as prefix for our script files, $galaxy is set in the startup script of the galaxy server
export MOVIEMAKERPATH="$galaxy""tools/customTools/movie_maker/"

#include library path to pymol in pythonpath, so python knows about the pymol 1.8.4 module
#include current directory in pythonpath, so scripts are available to import
export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${MOVIEMAKERPATH}:${PYTHONPATH}"
#used in pymolscript to save a seccond output
export POLAR_INTERACTION_FILENAME="polar_interaction_partners.txt"

#echo "working on $MOVIEMAKERPATH"
#echo "executing $MOVIEMAKERPATH""movie_maker_basic.py"
#echo "got "$#" arguments" > /home/webservices/philipp/movie_maker_basic.log
#check number of passed arguments, if we have 11, we have no cofactor, if 13 cofactor and color_carbon_cofactor
if [[ $# -eq 11 ]]
    then
        /home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker_basic.py" --input "$1" --ligand_name $2 --chain_name $3 --color_blind_friendly $4 --binding_site_radius $7 --check_halogen_interaction $8 --water_in_binding_site $9 --color_carbon "${10}" --session_export_version ${11} > /home/webservices/philipp/movie_maker_basic.log
    else
        /home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker_basic.py" --input "$1" --ligand_name $2 --chain_name $3 --color_blind_friendly $4 --binding_site_radius $7 --check_halogen_interaction $8 --water_in_binding_site $9 --color_carbon "${10}" --session_export_version ${11} --cofactor_name ${12} --color_carbon_cofactor ${13} > /home/webservices/philipp/movie_maker_basic.log
fi

# move created pymol session from current directory to output directory
mv basic_movie.pse $5
# move textfile from pymol-script to destination directory
mv $POLAR_INTERACTION_FILENAME $6
