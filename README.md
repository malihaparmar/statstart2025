[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Protein contact maps 

This project focuses on extracting structural information from proteins 
by computing distance matrices and identifying atomic contacts. 
We use PDB files as inputs, and generate contact maps and graph 
representations to visualize interactions between the residues.

## How does it work? 

1) Input: The project takes 3D protein structures from PDB files.
2) 2. Distance Matrix Calculation:
    * For each protein, it extracts the anchor atom which in this case is 
    C-alpha (CA) atoms.
    * It computes a Euclidean distance matrix 
    between all residues using their 3D coordinates.


## Acknowledgements 
Collaborator [@asquiciarino] (https://github.com/asquiciarino)
Mentor 
This is part of [statstart2025] (https://hsph.harvard.edu/fellowship-special-program/statstart/)