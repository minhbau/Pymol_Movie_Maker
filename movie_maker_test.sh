#!/bin/bash
#echo $1 # path to pdb file (input)
#echo $2 # ligand name
#echo $3 # chain
#echo $4 # path?

export PYTHONPATH="/home/webservices/philipp/special_pymol/modules:${PYTHONPATH}"
export MOVIEMAKERPATH="$galaxy""tools/customTools/movie_maker/"

# ser environment variable to enable rendering on the server with movie.produce
export FREEMOL=~/dev/freemol-trunk/freemol

arr=$(echo $4 | sed -e 's/\(.*\)\/.*/\1\//')

#path="/home/judith/galaxy-dist/tools/customTools"
#path+="/xbscore_output_pdf.py"
echo ${arr}

#echo "working on $MOVIEMAKERPATH"
#echo "executing $MOVIEMAKERPATH""movie_maker_basic.py"
/home/webservices/philipp/special_pymol/pymol -c -u $MOVIEMAKERPATH"movie_maker_basic.py" -- $1 $2 $2 $3 $4 $5 $galaxy > /home/webservices/philipp/movie_maker_basic.log

#cd ${arr}
#
#echo "test working?" > test.pse
#mv blubb.pse $5
mv basic_movie.pse $4

