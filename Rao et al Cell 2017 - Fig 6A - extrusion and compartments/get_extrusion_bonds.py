import sys
import os
import subprocess as sp
from numpy.random import random_integers, rand
import argparse
import numpy as np

# Bonds become fixed with probability according to in.probs.
# Loops reset upon hitting each other, but fixed bonds remain.
# Fixed bonds reset with probability resetprob.

def get_extrusion_bonds(probs='in.probs', steps=1000, nbonds=10, resetprob=0.0):
    print "Generating extrusion bond list..."

    # load sliding probabilities. If two columns, CTCF is oriented.
    f = open(probs, 'r')
    fwd_probs, rev_probs = [], []
    for line in f:
        tokens = map(float, line.split())
        assert len(tokens) == 1 or len(tokens) == 2
        if len(tokens) == 2:
            fwd_probs.append(tokens[0])
            rev_probs.append(tokens[1])
        else:
            fwd_probs.append(tokens[0])
            rev_probs.append(tokens[0])
    f.close()

    polylen = len(fwd_probs)

    # list of fixed anchors - locations on the polymer
    # list of fixed anchors, bool for each anchor
    fixed = np.zeros((nbonds, 2)).astype("bool")

    # list of bonds: timesteps x num_bonds x 2
    bonds = np.zeros((steps + 1, nbonds, 2)).astype("int")

    #### Initialize the bonds! ####
    # Initialize length=2 bonds at random positions, independent of boundaries.
    # Bond positions run from 1 to polylen
    start_pos = random_integers(3, polylen-2, nbonds)
    for i in range(nbonds):
        while np.sum(np.abs(start_pos - start_pos[i]) <= 2) >= 2:
            start_pos[i] = random_integers(3, polylen-2)
        bonds[0, i, 0] = start_pos[i] - 1
        bonds[0, i, 1] = start_pos[i] + 1


    #### slide the bonds! ####
    for step in range(steps):
        if (step+1) % 1000 == 0:
            print "step %i..." % (step+1)
        
        # shift all the extruding bonds
        for ndx in range(nbonds):
            bond = bonds[step, ndx]

            # reset bond?
            reset_bond = False
            # fix bonds?
            if rand() > fwd_probs[bond[0]-1]:
                fixed[ndx, 0] = True
            if rand() > rev_probs[bond[1]-1]:
                fixed[ndx, 1] = True

            # slide bonds
            nextbond = [max(bond[0] - 1, 2), \
                        min(bond[1] + 1, polylen-2)]

            # check for fixed loop anchors
            if fixed[ndx, 0]:
                nextbond[0] = bond[0]
            if fixed[ndx, 1]:
                nextbond[1] = bond[1]
            
            # RESET BOND?
            # list other bond anchors
            other_anchors = np.hstack([bonds[step+1, :ndx].flatten(),
                                       bonds[step, ndx+1:].flatten()])
            # check if nextbond conflicts with existing bonds
            if np.any(np.abs(other_anchors - nextbond[0]) == 0):
                reset_bond = True
            if np.any(np.abs(other_anchors - nextbond[1]) == 0):
                reset_bond = True
            # check endpoints
            if nextbond[0] == 2 or nextbond[1] == polylen-2:
                reset_bond = True
            
            # reset bond randomly?
            if rand() < resetprob:
                reset_bond = True

            if reset_bond:
                # unfix anchors
                fixed[ndx, :] = False
                retry = True
                while retry:
                    start_pos = random_integers(3, polylen-2)
                    retry = np.any(np.abs(other_anchors - start_pos) <= 1)
                nextbond[0] = start_pos - 1
                nextbond[1] = start_pos + 1

            # record shifted bond
            bonds[step + 1, ndx, 0] = nextbond[0]
            bonds[step + 1, ndx, 1] = nextbond[1]
        
    return bonds

def print_bonds(bonds, outfile='bond_history.txt'):
    steps, nbonds, _ = bonds.shape
    f = open(outfile, 'w')

    for step in range(steps-1):
        f.write("----- Sliding bonds -----\n")

        for b in range(nbonds):
            f.write("Shift:\t(%i, %i) --> (%i, %i)\n" % \
                    tuple(bonds[step, b].tolist() + bonds[step+1, b].tolist()))

    f.close()

        
