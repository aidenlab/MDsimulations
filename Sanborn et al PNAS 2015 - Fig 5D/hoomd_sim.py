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

############
# Initialize

# generate or load initial structure
if input_file:
    print "*** Loading initial coordinates from", input_file
    coords = np.loadtxt(input_file)
    assert len(coords) == int(params['N'])
else:
    # generate self-avoiding walk for initial conditions
    coords = get_poly_coords(int(params['N']))

# initialize system and box size
boxsize = np.max(np.abs(coords)) * 2.3
print "Boxsize:", boxsize
snapshot = data.make_snapshot(N=int(params['N']),
                              box=data.boxdim(L=boxsize),
                              particle_types=['A'],
                              bond_types=['bb', 'loop'])
system = init.read_snapshot(snapshot)

# initialize bead positions
ndx = 0
for ndx in range(len(system.particles)):
    p = system.particles[ndx]
    p.position = tuple(coords[ndx])
    
    if ndx > 0:
        p0 = system.particles[ndx-1]
        system.bonds.add("bb", p0.tag, p.tag)

# initialize bonds and forces
harmonic = bond.harmonic()
harmonic.bond_coeff.set('bb', k=params['k_bb'], r0=params['r0_bb'])
harmonic.bond_coeff.set('loop', k=params['k_loop'], r0=params['r0_loop'])

lj = pair.lj(r_cut=params['lj_cut'])
lj.pair_coeff.set('A', 'A', epsilon=params['epsilon'], sigma=params['sigma'])

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
if not input_file:
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


