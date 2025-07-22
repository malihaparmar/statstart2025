install.packages("bio3d")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("igraph")
install.packages("viridis")
library(bio3d)
library(reshape2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(igraph)
library(viridis)
my_pdb <- read.pdb("~/Downloads/6obi.pdb")
my_pdb2 <- read.pdb("~/Downloads/2jst.pdb")
my_pdb3 <- read.pdb("~/Downloads/3lel.pdb")

my_pdb$atom
my_pdb$seqrs
my_pdb$xyz

#summary : Euclidean distance between two 3D coordinate vectors
#Inputs: 
# v1, v2 - numeric vectors of length 3 (x, y, z coordinates)
#returns :  Euclidean distance between v1 and v2

euclidean_distance <- function(v1, v2) {
  if (length(v1) != length(v2)) {
    stop("Vector v1 and v2 must have the same length")
  }
  return(sqrt(sum((v1 - v2)^2)))
}

#summary : Computes a distance matrix using a specified atom 
#Inputs: 
#pdb - path to pdb file 
#anchor - atom type to use for distance calculation
#returns : matrix of Euclidean distances between anchor atoms

calc_dist_mat <- function(file, anchor = "CA") {
  pdb <- read.pdb(file)
  selected_idx <- which(pdb$atom$elety == anchor)
  
  filtered_table <- pdb$atom[selected_idx, ]
  if (nrow(filtered_table) == 0) stop("No atoms with specified anchor found.")
  
  coords <- filtered_table[, c("x", "y", "z")]
  n <- nrow(coords)
  dist_mat <- matrix(nrow = n, ncol = n)
  
  for (i in 1:n) {
    a1 <- coords[i, ]
    for (j in i:n) {
      a2 <- coords[j, ]
      dist <- euclidean_distance(a1, a2)
      dist_mat[i, j] <- dist
      dist_mat[j, i] <- dist
    }
  }
  return(dist_mat)
}

d1 <- calc_dist_mat("~/Downloads/6obi.pdb")
d2 <-  calc_dist_mat("~/Downloads/2jst.pdb")
d3 <- calc_dist_mat("~/Downloads/1e6v.pdb")

dim(d1)
dim(d2)

#summary :  Converts a distance matrix into a binary contact map 
#Inputs: 
#dist_matrix - distance matrix 
#threshold - distance threshold which residues are considered "in contact"
#returns : calculated binary matrix contact map

make_contact_map <- function (dist_matrix, threshold = 8.0) {
  dist_vec <- as.integer(dist_matrix < threshold) 
  return (matrix (dist_vec, nrow = nrow(dist_matrix), ncol = ncol(dist_matrix)))
}

make_contact_map(d1)
make_contact_map(d2)
make_contact_map(d3)

z <- make_contact_map(d1)
z2 <- make_contact_map(d2)
z3 <- make_contact_map(d3)

#summary : Converts a contact map to a graph object  
#Inputs: 
#contact map: binary matrix indicating residue-residue contacts 

contact_map_to_graph <- function (contact_map){
 adjacency_matrix <- graph_from_adjacency_matrix(contact_map, mode = "undirected")
 return(adjacency_matrix)
}

contact_map_to_graph(z)
contact_map_to_graph(z2)
contact_map_to_graph(z3)

#summary :Plots a graph representation of the contact network
#Inputs: 
# graph - igraph object created from contact map

plot_graph <- function (graph){
  plot(graph, layout = layout_with_kk(graph))
  
}

#summary : Plots a binary contact map as a heatmap
#Inputs: 
# dist_matrix - numeric matrix of residue-residue distances 


plot_dist_mat <- function(dist_matrix) {
  melt(dist_matrix)
  m <- melt(dist_matrix)
  ggplot(data = m, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(256))  
}

#summary: 
#Inputs:  
# contact_map - binary matrix of residue-residue contacts

plot_contact_map <- function(contact_map){
  m <- melt(contact_map)
  ggplot(data = m, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = "cividis", name = "Contact Strength") +
    #scale_fill_gradientn(colors = heat.colors(256)) 
    theme_minimal() +
    coord_fixed()
}

plot_contact_mat(z)
plot_contact_mat(z2)
plot_contact_mat(z3)
plot_dist_mat(d1)
plot_dist_mat(d2)
plot_dist_mat(d3)
plot_graph(contact_map_to_graph(z))
plot_graph(contact_map_to_graph(z2))
plot_graph(contact_map_to_graph(z3))
plot_contact_map(z)
plot_contact_map(z2)
plot_contact_map(z3)









