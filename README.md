# SanbornRaoPNAS2015
Scripts and molecular dynamics simulation code from the Aiden Lab. Developed by Adrian Sanborn.

The molecular dynamics code here will reproduce the simulations presented in various publications. The HOOMD-Blue package is required, available from http://glotzerlab.engin.umich.edu/hoomd-blue/. Note that HOOMD was rewritten in version 2.x and the code here is not compatible with these newer versions. This code is compatible with v1.3.3, which can be compiled from source.

Simulations can be run via the command

> hoomd hoomd_sim.py

hoomd_sim.py executes the main simulation commands. hoomd_utils.py contains some helper code. get_extrusion_bonds.py pre-computes the positions of each CTCF-cohesin-mediated extrusion bond at each time point in the simulation.

in.probs contains the forward and reverse probability, for each monomer, that the extrusion complex will slide past it (i.e. 1 - halting probability). It should always have the same number of lines as the polymer length N.

parameters.txt contains the major parameters, including:
  N - number of monomers, each monomer usually represents 1kb of chromatin
  steps - number of simulation steps to run during the initial collapse
  print_every - number of time steps between saving simulation snapshots
  nextrude - number of extrusion steps to run
  steps_per_extrude - number of simulation time steps for each extrusion step
  nbonds - number of extrusion bonds, randomly initialized
  T - temperature
  gamma - solvent viscosity
  epsilon, sigma, and lj_cut - describe the strength and range of the Lennard-Jones forces between polymer beads.
  k_bb and r0_bb - the strength and length of the polymer backbone bonds
  k_loop and r0_loop - the strength and length of the extrusion loop bonds
  dt - size of each time step 
Many parameters above are described in more detail in the HOOMD documentation.

For simulations with compartmentalization, compartments.txt contains the compartment designations as a list of the number of consecutive monomers (first column) of a designed type (second column). parameters.txt then includes the following:
  comp_spread - the variability of compartment boundary in different replicates, measured in the number of beads
  epsilon_intra - strength of LJ forces between beads of the same type
  epsilon-inter - strength of LJ forces between beads of different types

Contacts are counted as any pair of beads within a distance of 1.5. VMD can be used to visualize simulations by first important sim_topology.xml and then sim_trajectory.dcd.
