from argparse import ArgumentError

from pymol import cmd
from sys import argv
import os
from string import ascii_uppercase

#PATH TO CURRENT DIRECTORY
MOVIE_MAKER_PATH = os.environ.get('MOVIEMAKERPATH')

#SESSION_VERSION = round(float(argv[4]),2)
SESSION_VERSION = 1.2
SESSION_NAME = "basic_movie.pse"
NUMBER_OF_ARGUMENTS = 7

cmd.reinitialize()


#set session_export to be of desired version 
cmd.set("pse_export_version", SESSION_VERSION)


# load pdb file (first argument)
#renaming the *.dat file to *.pdb, loading with pymol and renaming for galaxy
os.rename(argv[1], argv[1].split(".")[0]+".pdb")
cmd.load(argv[1].split(".")[0]+".pdb")
os.rename(argv[1].split(".")[0]+".pdb", argv[1])


# Parse commandline arguments

PDB_NAME = argv[1].split(".")[0]
PDB_FILENAME = PDB_NAME + ".pdb"

print("loading %s"% (PDB_FILENAME,))

# movie_fade available at https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/movie_fade.py
cmd.do("run %sfade_movie.py"% (MOVIE_MAKER_PATH, ))

valid_amino_acid_3letter_codes = set("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR WAT SUL HEM".split(" "))

def parse_commandline_options():
    options = {}

    # passed filename is first passed argument, but flags are also in argv
    print("python file loaded")
    for i, ele in enumerate(argv):
        print("argv[%s]" % i, ele)

    #command line parameters are pdb-name, ligand_name, chain_name, output_filepath
    # session_version?

    # if len(argv) < NUMBER_OF_ARGUMENTS:
    #     raise ArgumentError("Not enough arguments supplied, only got %s, expected %s" % (len(argv), NUMBER_OF_ARGUMENTS))

    # resi = argv[3]
    ligand_name = argv[3]
    chain_name = argv[4]


    # Check whether we have a valid amino acid 3 letter code
    # if str(resi).upper() not in valid_amino_acid_3letter_codes:
    #     raise ArgumentError("Not a valid three letter amino acid code supplied.")

    #TODO create standard case in which no ligand name is required
    options["ligand_name"] = ligand_name
    options["chain_name"] = chain_name
    # options["path"] = path
    #output directory is managed by shellscript
    return options


def apply_colorblind_save(color_dict):
    #import color settings
    cmd.do("run %scolorblindfriendly.py" % MOVIE_MAKER_PATH)

    color_dict['protein_surface'] = "cb_sky_blue"
    color_dict['protein_cartoon'] = "cb_blue"
    color_dict['ligand'] = "cb_yellow"
    color_dict['cofactor'] = "cb_orange"
    color_dict['binding_site'] = "grey50"
    color_dict['interaction_polar'] = "cb_blue"
    # color_dict['oxygen'] = "cb_red"
    color_dict['oxygen'] = "cb_reddish_purple"
    color_dict['nitrogen'] = "cb_blue"

    return color_dict

def apply_settings(cmd_options):
    
    settings_dict = {}
    # BG Color is white
    cmd.bg_color('white')
    
    # Colors
    color_dict = apply_colorblind_save({})
    # color_dict['protein_surface'] = "lightblue"
    # color_dict['protein_cartoon'] = "skyblue"
    # color_dict['ligand'] = "yellow"
    # color_dict['cofactor'] = "sand"
    # color_dict['binding_site'] = "grey50"


    settings_dict["colors"] = color_dict
    settings_dict["cartoon_transparency"] = 0.6
    settings_dict["binding_site_radius"] = 5.0

    settings_dict["chain_name"] = cmd_options['chain_name']
    settings_dict["ligand_name"] = cmd_options['ligand_name']
    # settings_dict["path"] = cmd_options['path']

    """
    color_protein_surface = "lightblue"
    color_protein_cartoon = "skyblue"
    color_ligand = "yellow"
    color_cofactor = "sand"
    color_binding_site = "grey50"

    # Other settings
    cartoon_transparency = 0.6
    bs_radius = 5.0
    """

    return settings_dict


def create_selections(options):

    # Hide everything
    cmd.hide("lines")
    cmd.hide("nonbonded")

    # Protein structure
    cmd.select("protein_structure", "all")

    # Ligand
    #cmd.select("sele_ligand", "resn SUV")	#################################
    number_of_ligands_selected = cmd.select("sele_ligand", "organic and chain %s and resn %s" % (options['chain_name'], options["ligand_name"]))

    #Feature: If we did not get a correct chain name from the user, we will try to guess through the whole alphabet to find it
    # otherwise the ligand will not appear in the visualization
    if not number_of_ligands_selected:
        #list holding all ascii-letters
        for letter in ascii_uppercase:
            number_of_ligands_selected = cmd.select("sele_ligand", "organic and chain %s and resn %s" % (letter, options["ligand_name"]))
            if number_of_ligands_selected:
                break

    cmd.create("ligand", "sele_ligand")
    cmd.show("sticks", "ligand")
    cmd.color(options["colors"]['ligand'], "ligand and e. C")
    cmd.color(options["colors"]['oxygen'], "ligand and e. O")
    cmd.color(options["colors"]['nitrogen'], "ligand and e. N")

    # Cofactor
    #cmd.select("sele_cofactor", "resn OLA")	#################################
    #cmd.create("cofactor", "sele_cofactor")
    #cmd.show("sticks", "cofactor")
    #cmd.color(color_cofactor, "cofactor and e. C")

    # Surface
    cmd.create("protein_surface", "all")
    cmd.hide("lines", "protein_surface")
    cmd.hide("sticks", "protein_surface")
    cmd.hide("nonbonded", "protein_surface")
    cmd.show("surface", "protein_surface")
    cmd.color(options["colors"]['protein_surface'], "protein_surface")

    # Cartoon
    cmd.copy("protein_cartoon", "protein_surface")
    cmd.hide("surface", "protein_cartoon")
    cmd.show("cartoon", "protein_cartoon")
    cmd.color(options["colors"]['protein_cartoon'], "protein_cartoon")

    # Transparent Cartoon
    cmd.copy("protein_cartoon_transparent", "protein_cartoon")
    cmd.set("cartoon_transparency", options['cartoon_transparency'], "protein_cartoon_transparent")

    # 2nd Transparent Cartoon
    cmd.copy("protein_cartoon_transparent_more", "protein_cartoon")
    cmd.set("cartoon_transparency", 0.9, "protein_cartoon_transparent_more")


    # Binding Site
    cmd.select("sele_binding_site", "ligand expand %s"%(options['binding_site_radius']))
    cmd.select("sele_binding_site", "br. sele_binding_site")
    cmd.select("sele_binding_site", "sele_binding_site and protein_structure")
    cmd.select("sele_binding_site", "sele_binding_site and not sele_ligand")
    cmd.create("binding_site", "sele_binding_site")
    cmd.hide("surface", "binding_site")
    cmd.hide("nonbonded")
    cmd.show("sticks", "binding_site")
    cmd.show("nb_spheres", "binding_site")
    cmd.color(options["colors"]['binding_site'], "binding_site and e. C")
    cmd.color(options["colors"]['nitrogen'], "binding_site and e. N")
    cmd.color(options["colors"]['oxygen'], "binding_site and e. O")

    # TODO quickfix for basic movie
    cmd.select("interacting_residues", "sele_binding_site")

    # resi_list = [str(tup[1]) for tup in resn_tuples]
    # print resi_list
    # resi_list = [110, 111, 227, 324, 350, 4025]
    # resi_list = [str(asd) for asd in resi_list]
    # s = " or resi "
    # resi_or_string = s.join(resi_list)
    # non_resi_string = " and not resi ".join(resi_list)
    # print(resi_or_string)
    # print non_resi_string

    # cmd.create("interacting_residues", "sele_binding_site and resi %s" % (resi_or_string))
    # cmd.create("non_interacting_residues", "sele_binding_site and not resi %s" % (non_resi_string))
    # cmd.show("sticks", 'non_interacting_residues')
    # cmd.hide('non_interacting_residues')
    # cmd.hide('binding_site')
    # cmd.color('lightblue','interacting_residues')
    cmd.show("sticks", "interacting_residues")
    cmd.color("grey50", "interacting_residues and e. C")
    cmd.color(options["colors"]['nitrogen'], "interacting_residues and e. N")
    cmd.color(options["colors"]['oxygen'], "interacting_residues and e. O")

    # Polar contacts
    cmd.distance("interaction_polar", "ligand", "binding_site", mode=2)
    cmd.hide("labels", "interaction_polar")
    cmd.color(options['colors']["interaction_polar"], "interaction_polar")

    # other polar interactions
    # cmd.distance("interaction_polar", "sele_interacting_HIS350", "sele_interacting_HOH4025", mode=2)
    # additional hbb to water in pocket


    # Halogen Bond
    # cmd.select("sele_chlorine", "e. Cl and ligand")
    # cmd.distance("interaction_halogen_bond_distance", "sele_chlorine", "sele_interacting_ALA110 and n. O")	############
    # cmd.angle("interaction_halogen_bond_angle", "neighbor sele_chlorine", "sele_chlorine", "sele_interacting_ALA110 and n. O")
    # cmd.color("gold", "interaction_halogen_bond_distance")
    # cmd.color("gold", "interaction_halogen_bond_angle")
    # cmd.hide("labels")


    # Add more interactions using names such as "interaction_halogen_bond" or "inteteraction_CH_pi"


def do_it():
    # run all script components
    commandline_options = parse_commandline_options()
    settings_dict = apply_settings(commandline_options)
    create_selections(settings_dict)
    create_views(settings_dict)
    #create scenes and frames for movie, execute a pymol script with @
    cmd.do("@%smovie_maker_basic_script.pml" % MOVIE_MAKER_PATH)
    #Save session
    cmd.save("basic_movie.pse")



def print_binding_site_residues():
    # print all binding site residues with residue numbers
    model_binding_site = cmd.get_model("binding_site")
    my_resns_set = set([])
    my_resns = []
    for resid in model_binding_site.atom:
        a_tuple = (resid.resn, resid.resi)
        resn, resi = a_tuple
        if a_tuple not in my_resns_set:
            my_resns_set.add(a_tuple)
            #print a_tuple
            #print('cmd.select("sele_interacting_%s%s", "sele_binding_site and resi %s")' % (resn,resi,resi))
            my_resns.append(a_tuple)
    return my_resns

#resn_tuples = print_binding_site_residues()

#cmd.select("sele_interacting_ALA110", "sele_binding_site and resi 110")
#cmd.select("sele_interacting_THR111", "sele_binding_site and resi 111")
#cmd.select("sele_interacting_PHE227", "sele_binding_site and resi 227")
#cmd.select("sele_interacting_ASN324", "sele_binding_site and resi 324")
#cmd.select("sele_interacting_HIS350", "sele_binding_site and resi 350")
#cmd.select("sele_interacting_HOH4025", "sele_binding_site and resi 4025")


def create_views(options):

    cmd.disable("all")
    cmd.orient("protein_surface")
    #turn 90 degrees around x and z
    #orient correctly, top accessible region
    #cmd.turn("z", -90)
    cmd.zoom("protein_surface", 5)
    cmd.enable("protein_surface")
    cmd.set("transparency", 0.5)
    cmd.enable("protein_cartoon")
    cmd.view("1", action="store")
    cmd.scene("F1", action="store")


    # 2 and F2 show ligand in pocket
    cmd.disable("all")
    cmd.enable("protein_surface")
    cmd.set("transparency", 0.5)
    cmd.enable("protein_cartoon")
    cmd.enable("ligand")
    cmd.view("2", action="store")
    cmd.scene("F2", action="store")


    # 3 and F3
    cmd.disable("all")
    cmd.enable("protein_cartoon")
    cmd.enable("ligand")
    cmd.view("3", action="store")
    cmd.scene("F3", action="store")


    # 4 and F4
    cmd.disable("all")
    cmd.enable("ligand")
    cmd.view("4", action="store")
    cmd.scene("F4", action="store")

    # 5 and F5
    cmd.disable("all")
    cmd.enable("ligand")
    cmd.enable("protein_cartoon")
    cmd.orient("ligand")
    cmd.zoom("interacting_residues", 2)
    cmd.view("5", action="store")
    cmd.scene("F5", action="store")

    # 6 and F6
    cmd.disable("all")
    cmd.enable("ligand")
    cmd.enable("interacting_residues")
    cmd.enable("interaction_polar")
    cmd.zoom("interacting_residues", 2)
    #cmd.set_view("""\
    #                 -0.652670681,   -0.733961642,    0.187926084,\
    #                  0.654971540,   -0.671278477,   -0.346988708,\
    #                  0.380826712,   -0.103380904,    0.918851793,\
    #                 -0.000000000,    0.000000000,  -73.133735657,\
    #                 52.917945862,    8.675487518,   53.204895020,\
    #                 55.469581604,   90.797828674,  -20.000000000 """)
    cmd.view("6", action="store")
    cmd.scene("F6", action="store")

    # 7 and F7
    cmd.disable("all")
    cmd.enable("ligand")
    cmd.enable("binding_site")
    cmd.enable("interacting_residues")
    cmd.enable("interaction_polar")
    cmd.zoom("interaction_polar", 2)
    #cmd.enable("interaction_halogen_bond_distance")
    #cmd.enable("interaction_halogen_bond_angle")
    #cmd.hide("sticks", "binding_site")
    #cmd.hide("nb_spheres", "binding_site")
    #cmd.set_view ("""\
    #        -0.685239613,   -0.717070580,   -0.127501234,\
    #         0.620792449,   -0.483501613,   -0.617124677,\
    #         0.380877644,   -0.502027392,    0.776465952,\
    #         0.000000000,   -0.000000000,  -98.679916382,\
    #        51.064998627,    6.873000145,   55.061000824,\
    #        84.767715454,  112.592117310,  -20.000000000 """)
    cmd.view("7", action="store")
    cmd.scene("F7", action="store")
    cmd.viewport(500, 500)

'''
#record the movie
#run script with @script.txt
#cmd.viewport(2000, 2000)
"""
mset 1x2000
mview store, 1, scene=F1
turn y, 120
mview store, 100, power = 1.0
turn y, 120[
mview store, 200, power = 1.0
mview store, 300, scene=F1
turn x, 120
mview store, 400, power = 1.0
turn x, 120
mview store, 500, power = 1.0
mview store, 600, scene=F1
mview store, 601, scene=F2
turn x, 100
mview store, 700, power = 1.0
turn x, -100
mview store, 800, power = 1.0
mview store, 900, scene=F3
turn y, 120
mview store, 1000, power = 1.0
turn y, 120
mview store, 1100, power = 1.0
mview store, 1200, scene=F3
movie_fade cartoon_transparency, 1201, 0.0, 1290, 1.0
mview store, 1300, scene=F5
turn x, 50
mview store, 1400, power = 1.0
turn x, 50
mview store, 1500, scene=F5
mview store, 1600, scene=F6
mview store, 1700, scene=F7
turn y, 60
mview store, 1750, power = 1.0
turn y, -120
mview store, 1850, power = 1.0
mview store, 1900, scene=F7
mview store, 2000, scene=F7
set ray_trace_frames, 1
"""


"""
#TODO get all polar interacting amino acids and zoom into each of them
# 10 frames per AA
mset 1 x1440
mview store

# this code maps the zooming of
# one AA and it's neighbor to 10 frames
python
for x in range(0,144):
  cmd.frame((10*x)+1)
  cmd.zoom( "n. CA and i. " + str(x) + "+" + str(x+1))
  cmd.mview("store")
python end
"""
'''

do_it()

