[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Protein contact maps 

This project focuses on extracting structural information from proteins 
by computing distance matrices and identifying atomic contacts. 
We use PDB files as inputs, and generate contact maps and graph 
representations to visualize interactions between the residues.

## How does it work? 

1) The project extracts 3D protein structures specifically a residue's x, y, z, 
from PDB files.
    * Function used: read.pdb() from the bio3d package
2) Distance Matrix Calculation:
    * Function used: calc_dist_mat(file, anchor = "CA")
    * Internally uses: euclidean_distance(v1, v2) to calculate 3D distances
    * For each protein, it extracts the anchor atom which in this case is 
    C-alpha (CA) atoms.
    * A Euclidean distance matrix is computed for all pairs of residues using
    their 3D coordinates.
    * To optimize performance, only the upper triangle of the matrix is 
    computed, and the values are mirrored to the lower triangle (since the 
    matrix is symmetric).
3) Contact Map: 
  * The distance matrix is converted into a binary contact map using a cutoff 
  threshold of 8 Ångströms.
  * If the distance between two residues is ≤ 8Å, they are marked as being
  "in contact" (1); otherwise, not in contact (0).
  * This contact map provides a simplified view of spatial residue interactions.
4) Graphs: 
  * Heatmap - heat map is graphed of the pairwise distance matrix  
  * Histogram - plotted for the degree distribution of protein conatcts  
  * Network graph - plotted for the amino acid contacts 


## Acknowledgements 
Collaborator [@asquiciarino] (https://github.com/asquiciarino)
Mentor 
This is part of [statstart2025] (https://hsph.harvard.edu/fellowship-special-program/statstart/)