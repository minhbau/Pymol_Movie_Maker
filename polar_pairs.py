'''
Authors:
    2011 Thomas Holder, MPI for Developmental Biology
    Available at https://pymolwiki.org/index.php/Polarpairs

    2017 Philipp Weber, philmaweb@gmail.com


Example usage:

pairs = polarpairs('chain A', 'chain B')
for p in pairs:
    dist = cmd.get_distance('(%s`%s)' % p[0], '(%s`%s)' % p[1])
    print p, 'Distance: %.2f' % (dist)

'''

from pymol import cmd

def polarpairs(sel1, sel2, cutoff=4.0, angle=63.0, name='', state=1, quiet=1):
    '''
ARGUMENTS

    sel1, sel2 = string: atom selections

    cutoff = float: distance cutoff

    angle = float: h-bond angle cutoff in degrees. If angle="default", take
    "h_bond_max_angle" setting. If angle=0, do not detect h-bonding.

    name = string: If given, also create a distance object for visual representation

SEE ALSO

    cmd.find_pairs, cmd.distance
    '''
    cutoff = float(cutoff)
    quiet = int(quiet)
    state = int(state)
    if angle == 'default':
        angle = cmd.get('h_bond_max_angle', cmd.get_object_list(sel1)[0])
    angle = float(angle)
    mode = 1 if angle > 0 else 0
    x = cmd.find_pairs('(%s) and donors' % sel1, '(%s) and acceptors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle) + \
        cmd.find_pairs('(%s) and acceptors' % sel1, '(%s) and donors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle)
    x = sorted(set(x))
    if not quiet:
        print 'Settings: cutoff=%.1fangstrom angle=%.1fdegree' % (cutoff, angle)
        print 'Found %d polar contacts' % (len(x))
    if len(name) > 0:
        for p in x:
            cmd.distance(name, '(%s`%s)' % p[0], '(%s`%s)' % p[1])
    return x

cmd.extend('polarpairs', polarpairs)



def polartuples(the_polarpairs, residue_name='polar_interaction'):
    '''
    get list of polar interacting residues from polarpair
    '''
    polar_resn_tuples = []
    polar_resn_tuple_set = set([])
    #we are only interested in the residues in the binding site, polarpairs returns pairs
    pairs = list(set([p[0] for p in the_polarpairs]))  # remove duplicate entries
    for i, p in enumerate(pairs):
        #create the object for our visualization
        pair_name = residue_name + '_%s' % (i,)
        sele_pair_name = "sele_" + residue_name + '_%s' % (i,)
        cmd.select(sele_pair_name, "(%s`%s)" % p)
        cmd.select(sele_pair_name, 'byres ' + sele_pair_name)
        cmd.create(pair_name, sele_pair_name)
        cmd.delete(sele_pair_name)
        cmd.show('sticks', pair_name)
        temp_model = cmd.get_model(pair_name)

        #extract resn, resi, and chain identifier from binding site
        for atom in temp_model.atom:
            a_tuple = (atom.resi, atom.resn, atom.chain)
            #remove all duplicate entries
            if a_tuple not in polar_resn_tuple_set:
                polar_resn_tuples.append(a_tuple)
                polar_resn_tuple_set.add(a_tuple)
                print a_tuple
                break
        # print p, 'Distance: %.2f' % (dist)

    return polar_resn_tuples

cmd.extend('polartuples', polartuples)


