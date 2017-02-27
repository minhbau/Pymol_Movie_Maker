#!/bin/bash
#execute our lovely script with the correct pymol and passed commondaline parameters
#activate conda environment
source activate special_pymol
#$1 is path to pdb file used for the visualization
pymol -c -u $galaxy$"tools/customTools/movie_maker/movie_maker_basic.py" -- $1
