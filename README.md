# dijkstra-water-wire
This is a C / C++ program for Dijkstra's single source shortest
path algorithm. The program is for adjacency matrix
representation of the graph.
taken from http://www.geeksforgeeks.org/printing-paths-dijkstras-shortest-path-algorithm/ .
The code has been modified to calculate water wires connecting a donor and acceptor from an MD simulation.
Currently, only Amber MD has been tested and the code can use either a pdb trajectory file 
or an Amber ASCII trajectory file (mdcrd). For a path to be accepted, waters must be within 
a user defined cutoff of each other
other or of the donor or acceptor.
=====================
Anticipated code updates:
1) Allowing for angles. In other words, viability also has to do with the Donor-hydrogen-acceptor angle.
2) The code can do pdb *or* mdcrd. Make one code that can do both given user input.
============================
Anticipated file updates
1) Example pdb/mdcrd
2) example outputs
