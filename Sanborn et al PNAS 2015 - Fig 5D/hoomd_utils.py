#!/usr/bin/env python

from hoomd_script import *
import numpy as np

# Read parameter file in the simulation folder
def read_params(param_file='parameters.txt'):
    print "*** Loading parameters"
    f = open(param_file, 'r')
    params = {}

    for line in f:
        if '#' in line: continue # comments
        if '=' not in line: continue

        tokens = line.split('=')
        assert len(tokens) == 2
        params[tokens[0].strip()] = float(tokens[1].strip())
        print tokens[0].strip(), tokens[1].strip()

    return params


# Generates self-avoiding walk of length N using HOOMD functions
def get_poly_coords(N):
    done = False
    tries = 0 
    Lmax = N/10

    while not done and tries < 10:
        print "*** Generating polymer coordinates, try %i of 10" % tries

        # Generate polymer coordinates using HOOMD functions.
        # Use bond lengths 1/10th the size to save memory usage.
        polymer = dict(bond_len=0.1, type=['A']*N, bond="linear", count=1)
        system = init.create_random_polymers(box=data.boxdim(L=Lmax),
                                             polymers=[polymer],
                                             separation=dict(A=0.0499),
                                             seed=np.random.randint(1e5))

        # read out coordinates, shift each axis so that no wrapping occurs
        coords = [p.position for p in system.particles]
        coords = np.mod(np.array(coords), Lmax*2)
        for D in range(3):
            while np.min(coords[:, D]) < 10:
                shift = np.array([0, 0, 0])
                shift[D] = -Lmax/10
                coords = np.mod(coords + shift, Lmax*2)

        # rescale and shift coordinates
        coords *= 10
        coords -= np.mean(coords, axis=0)

        del system
        init.reset()

        # verify that polymer is not wrapped around the box boundaries
        done = np.max(np.abs(coords[:-1, :] - coords[1:, :])) < 2.0
        tries += 1

    if done:
        return coords
    else:
        raise NotImplementedError

