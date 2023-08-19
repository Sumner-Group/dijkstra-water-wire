# dijkstra-water-wire
This is a C / C++ program for Dijkstra's single source shortest
path algorithm based off of code from http://www.geeksforgeeks.org/printing-paths-dijkstras-shortest-path-algorithm/, by Aditya Goel.
The code has been modified to calculate water wires connecting a donor and acceptor from an MD simulation.
Currently, only Amber MD has been tested and the code can use either a pdb trajectory file 
or an Amber ASCII trajectory file (mdcrd). For a path to be accepted, waters must be within 
a user defined cutoff of each other
other or of the donor or acceptor. 

Input file description

filetype (mdcrd/pdb) 
filename 
number of atoms 
atom number for start of water chain 
atom number for end of water chain 
number of water molecules  (**Only for mdcrd files**)
distance cut-off

To compile the code

g++ -O2 -o path.exe water-wire.cpp

To run the code

./path.exe < inputfile

Sample inputs are given as input.pdb and input.mdcrd

Sample outputs are given as output.pdb and output.mdcrd

Sample mdcrd and pdb files are given as ub-nvt-step1-3.pdb and ub-nvt-step13.mdcrd. They each have three geometries.

=====================

Anticipated updates:
1) Allowing for angles. In other words, viability also has to do with the donor-hydrogen-acceptor angle.

**Do you want to cite this work?**
---
Wilson, R. H.; Zamfir, S.; Sumner, I. “Molecular dynamics simulations reveal a new role for a conserved active site asparagine in a ubiquitin-conjugating enzyme” Journal of Molecular Graphics and Modelling, 2017, 76, 403-411. [doi:10.1016/j.jmgm.2017.07.006](https://doi.org/10.1016/j.jmgm.2017.07.006)

[![DOI](https://zenodo.org/badge/86590973.svg)](https://zenodo.org/badge/10.5281/zenodo.495640/86590973)
