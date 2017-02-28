from pymol import cmd
from sys import argv
import os
#from shutil import copyfile

#command line parameters are pdb-name, chain_name, ligand_name, session_version, output_filepath
# passed filename is first passed argument, but flags are also in argv
#pdb_name = argv[-2].split(".")[0]
pdb_name = argv[1].split(".")[0]
#print "all arguments", argv
#print pdb_name
pdb_filename = pdb_name + ".pdb"

#output directory is managed by shellscript
#output_filename = argv[-1]
path_to_galaxy_dir = argv[6]

#original rename script from show_ligand
os.rename(argv[1], argv[1].split(".")[0]+".pdb")
cmd.load(argv[1].split(".")[0]+".pdb")
os.rename(argv[1].split(".")[0]+".pdb", argv[1])

#chain_name = argv[]
#ligand_name = argv[]
#ligand_name = argv[]
#session_version = round(float(argv[4]),2)
session_version = 1.2
#session_name = "basic_movie_%s_%s.pse" % (pdb_name, session_version)
session_name = "basic_movie_%s.pse" % (session_version,)

pathname = os.path.dirname(argv[0])
os.path.realpath(__file__)
full_path = os.path.abspath(pathname)+"/%s"% (pdb_filename,)
print('full path =', full_path)



cmd.reinitialize()
cmd.do("run fade_movie.py")
# movie_fade available at https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/movie_fade.py
cmd.bg_color("white")

#set session_export to be of desired version #TODO
cmd.set("pse_export_version", session_version)

#copy .dat file to .pdb
#os.rename(argv[-1], pdb_filename)
#copyfile(argv[-2], pdb_filename)
#load file into pymol
cmd.load(pdb_filename)

# Colors
color_protein_surface = "lightblue"
color_protein_cartoon = "skyblue"
color_ligand = "yellow"
color_cofactor = "sand"
color_binding_site = "grey50"

# Other settings
cartoon_transparency = 0.6
bs_radius = 5.0

# Hide everything
cmd.hide("lines")
cmd.hide("nonbonded")

# Protein structure
cmd.select("protein_structure", "all")

# Ligand
#"ligand"
#cmd.select("sele_ligand", "resn SUV")	#################################
cmd.select("sele_ligand", "organic")	#################################
cmd.create("ligand", "sele_ligand")
cmd.show("sticks", "ligand")
cmd.color(color_ligand, "ligand and e. C")

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
cmd.color(color_protein_surface, "protein_surface")

# Cartoon
cmd.copy("protein_cartoon", "protein_surface")
cmd.hide("surface", "protein_cartoon")
cmd.show("cartoon", "protein_cartoon")
cmd.color(color_protein_cartoon, "protein_cartoon")

# Transparent Cartoon
cmd.copy("protein_cartoon_transparent", "protein_cartoon")
cmd.set("cartoon_transparency", cartoon_transparency, "protein_cartoon_transparent")

# 2nd Transparent Cartoon
cmd.copy("protein_cartoon_transparent_more", "protein_cartoon")
cmd.set("cartoon_transparency", 0.9, "protein_cartoon_transparent_more")


# Binding Site
cmd.select("sele_binding_site", "ligand expand %s"%(bs_radius))
cmd.select("sele_binding_site", "br. sele_binding_site")
cmd.select("sele_binding_site", "sele_binding_site and protein_structure")
cmd.select("sele_binding_site", "sele_binding_site and not sele_ligand")
cmd.create("binding_site", "sele_binding_site")
cmd.hide("surface", "binding_site")
cmd.hide("nonbonded")
cmd.show("sticks", "binding_site")
cmd.show("nb_spheres", "binding_site")
cmd.color(color_binding_site, "binding_site and e. C")



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

#TODO quickfix for basic movie
cmd.select("interacting_residues", "sele_binding_site")

#resi_list = [str(tup[1]) for tup in resn_tuples]
#print resi_list
#resi_list = [110, 111, 227, 324, 350, 4025]
#resi_list = [str(asd) for asd in resi_list]
#s = " or resi "
#resi_or_string = s.join(resi_list)
#non_resi_string = " and not resi ".join(resi_list)
#print(resi_or_string)
#print non_resi_string

#cmd.create("interacting_residues", "sele_binding_site and resi %s" % (resi_or_string))
#cmd.create("non_interacting_residues", "sele_binding_site and not resi %s" % (non_resi_string))
#cmd.show("sticks", 'non_interacting_residues')
#cmd.hide('non_interacting_residues')
#cmd.hide('binding_site')
#cmd.color('lightblue','interacting_residues')
cmd.show("sticks", "interacting_residues")
cmd.color("grey50", "interacting_residues and e. C")


# Polar contacts
cmd.distance("interaction_polar", "ligand", "binding_site", mode=2)
cmd.hide("labels", "interaction_polar")
cmd.color("blue", "interaction_polar")

#other polar interactions
#cmd.distance("interaction_polar", "sele_interacting_HIS350", "sele_interacting_HOH4025", mode=2)
#additional hbb to water in pocket


# Halogen Bond
#cmd.select("sele_chlorine", "e. Cl and ligand")
#cmd.distance("interaction_halogen_bond_distance", "sele_chlorine", "sele_interacting_ALA110 and n. O")	############
#cmd.angle("interaction_halogen_bond_angle", "neighbor sele_chlorine", "sele_chlorine", "sele_interacting_ALA110 and n. O")
#cmd.color("gold", "interaction_halogen_bond_distance")
#cmd.color("gold", "interaction_halogen_bond_angle")
#cmd.hide("labels")


# Add more interactions using names such as "interaction_halogen_bond" or "inteteraction_CH_pi"

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

cmd.save(session_name)
#rock 40 degrees up and down

#record the movie
"""
cmd.viewport(2000, 2000)
mset 1x2000
mview store, 1, scene=F1
turn y, 120
mview store, 100, power = 1.0
turn y, 120
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
cmd.set("ray_trace_frames", 1)
"""




