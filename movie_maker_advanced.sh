#!/bin/bash
#execute our lovely script with the correct pymol, python library support and passed commandline parameters
#echo $1  # path to pdb file (input)
#echo $2  # ligand name
#echo $3  # chain
#echo $4  # use colorblind friendly coloring?
#echo $5  # path pse file
#echo $6  # path to polar interaction .txt file
#echo $7  # path pymol-script .pml file
#echo $8  # binding_site_radius
#echo $9  # check_halogen_interaction
#echo $10  # water_in_binding_site
#echo $11 # color_carbon
#echo $12 # session_export_version
#echo $13 # color of polar-interactions
#echo $14 # cofactor name
#echo $15 # color of carbon in cofactor

#used in our pymolscript as prefix for our script files, $galaxy is set in the startup script of the galaxy server
export MOVIEMAKERPATH="$galaxy""tools/customTools/movie_maker/"

#include library path to pymol in pythonpath, so python knows about the pymol 1.8.4 module
#include current directory in pythonpath, so scripts are available to import
export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${MOVIEMAKERPATH}:${PYTHONPATH}"
#used in pymolscript to save a seccond output
export POLAR_INTERACTION_FILENAME="polar_interaction_partners.txt"
export MOVIE_SCRIPT_FILENAME="${MOVIEMAKERPATH}""movie_maker_basic_script.pml"

#echo "working on $MOVIEMAKERPATH"
#echo "executing $MOVIEMAKERPATH""movie_maker_basic.py"
#echo "got "$#" arguments" > /home/webservices/philipp/movie_maker_basic.log
#check number of passed arguments, if we have 13, we have no cofactor, if 15 cofactor and color_carbon_cofactor
if [[ $# -eq 13 ]]
    then
        /home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker.py" --input "$1" --ligand_name $2 --chain_name $3 --color_blind_friendly $4 --binding_site_radius $8 --check_halogen_interaction $9 --water_in_binding_site "${10}" --color_carbon "${11}" --session_export_version ${12} --color_polar_interactions ${13} > /home/webservices/philipp/movie_maker_basic.log
    else
        /home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker.py" --input "$1" --ligand_name "$2" --chain_name "$3" --color_blind_friendly "$4" --binding_site_radius "$8" --check_halogen_interaction $9 --water_in_binding_site "${10}" --color_carbon "${11}" --session_export_version ${12} --color_polar_interactions ${13} --cofactor_name ${14} --color_carbon_cofactor ${15} > /home/webservices/philipp/movie_maker_basic.log
fi

# move created pymol session from current directory to output directory
mv basic_movie.pse "$5"

# move textfile from pymol-script to destination directory
mv $POLAR_INTERACTION_FILENAME "$6"

# move pml script to destination directory, in combination with session allows easy modification of the movie
mv $MOVIE_SCRIPT_FILENAME "$7"