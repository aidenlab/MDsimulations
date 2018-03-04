#!/usr/bin/env python

import numpy as np
from hoomd_script import *
import sys
from hoomd_utils import read_params, get_poly_coords
from get_extrusion_bonds import get_extrusion_bonds, print_bonds

if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = None

context.initialize()

params = read_params("parameters.txt")
N = int(params['N'])

############
# Initialize

# generate or load initial structure
if input_file:
    print "*** Loading initial coordinates from", input_file
    coords = np.loadtxt(input_file)
    assert len(coords) == N
else:
    # generate self-avoiding walk for initial conditions
    coords = get_poly_coords(N)

# read compartent file
f = open('compartments.txt', 'r')
comps = []
comptypes = []
# spread the exact position of the compartment shift
if 'comp_spread' in params:
    spread = int(params['comp_spread'])
else:
    spread = 0
prev_adjust = 0
for line in f:
    tokens = line.split()
    assert len(tokens) == 2
    numof, type = tokens

    if type not in comptypes:
        comptypes.append(type)

    adjust = np.random.randint(-1*spread, spread+1)
    comps += [type] * (int(numof) - prev_adjust + adjust)
    print "Compartment boundary at", len(comps)
    prev_adjust = adjust # add/subtract shift to the next compartment boundary

# adjust the final boundary
if len(comps) > N:
    comps = comps[:N]
if len(comps) < N:
    comps += [type] * (N - len(comps))
print "Adjusted final boundary to", len(comps)
assert len(comptypes) == 2 # only support two types right now
f.close()

# initialize system and box size
boxsize = np.max(np.abs(coords)) * 10
print "Boxsize:", boxsize
snapshot = data.make_snapshot(N=N,
                              box=data.boxdim(L=boxsize),
                              particle_types=comptypes,
                              bond_types=['bb', 'loop'])
system = init.read_snapshot(snapshot)

# initialize bead positions
ndx = 0
for ndx in range(len(system.particles)):
    p = system.particles[ndx]
    p.position = tuple(coords[ndx])
    p.type = comps[ndx]
    
    if ndx > 0:
        p0 = system.particles[ndx-1]
        system.bonds.add("bb", p0.tag, p.tag)

# initialize bonds
harmonic = bond.harmonic()
harmonic.bond_coeff.set('bb', k=params['k_bb'], r0=params['r0_bb'])
harmonic.bond_coeff.set('loop', k=params['k_loop'], r0=params['r0_loop'])

# initialize inter-monomeric forces
lj = pair.lj(r_cut=params['lj_cut'])
lj.pair_coeff.set(comptypes[0], comptypes[0],
                  epsilon=params['epsilon_intra'], sigma=params['sigma'])
lj.pair_coeff.set(comptypes[1], comptypes[1],
                  epsilon=params['epsilon_intra'], sigma=params['sigma'])
lj.pair_coeff.set(comptypes[0], comptypes[1],
                  epsilon=params['epsilon_inter'], sigma=params['sigma'])

# initialize dynamics
all = group.all();
integrate.mode_standard(dt=params['dt'])
lang = integrate.langevin(group=all, T=params['T'], seed=5)
lang.set_gamma('A', gamma=params['gamma'])
#integrate.nve(group=all)

# initialize file output
xml = dump.xml(filename='sim_topology.xml', vis=True) # system topology
dump.dcd(filename='sim_trajectory.dcd', period=params['print_every']) # trajectory


###############
# Run LJ collapse
#if not input_file:
#    run(params['steps']+1)
run(params['steps']+1)


###############
# Extrusion
extr_bonds = get_extrusion_bonds(probs = 'in.probs',
                                 steps = int(params['nextrude']),
                                 nbonds = int(params['nbonds']))
print_bonds(extr_bonds)

for sndx in range(int(params['nextrude'])):
    # delete existing loop bonds
    loops = [bnd for bnd in system.bonds if bnd.type=='loop']
    for bnd in loops:
        system.bonds.remove(bnd.tag)

    # set new loop bonds
    for bndx in range(int(params['nbonds'])):
        b1, b2 = extr_bonds[sndx, bndx]
        system.bonds.add("loop", b1-1, b2-1)

    run(params['steps_per_extrude'])


