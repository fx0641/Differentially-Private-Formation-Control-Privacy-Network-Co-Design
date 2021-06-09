# Differentially Private Formation Control: Privacy and Network Co-Design

This repository contains a MATLAB implementation of the privacy and network co-design problem presented in the publication Differentially Private Formation Control: Privacy and Network Co-Design that is under review at IEEE Transacions on Control of Network Systems. In short, the algorithm here solves an optimization problem to compute the best privacy parameters and edge weights in a communication topology for a network of agents running differentially private formation control. Below you will find a brief description of each of the files in this repositiory, instructions for intended usage, as well as some sample outputs.

### Usage


This code was developed with MATLAB 2019a, I am not sure exactly what versions it is compatabile with. The code is designed so everything can be adjusted and run with 'main.m'. The first section of this code is titled inputs, in this section you can provide numerical values for all of your parameters and specifiy an input graph Laplacian. To test the solver, I also created a function to generate random graph laplacians given a number of agents $N$. This function is called 'makeRandomGraph.m' and just takes $N$ as an input.

Since we are working with undirected graphs over $N$ nodes, we can exploit symmetry and just optimize over the elements in the upper triangular part of the graph laplacian and impose symmetry. This is used wherever possible in the code to maximize effeciency.

### Description

 - main.m: Initializes everything and calls the other functions
 - objectiveFunction.m: The objective function for the optimization problem
 - constraints.m: The constraints for the optimization problem
 - makeRandomGraph.m: Makes a random graph Laplacian over N nodes
 - makeLaplacian.m: Takes the upper triangular elements of a graph laplacian and makes a graph laplacian
 - calculateSSerror.m: Calculates the steady state error bound found in the paper.

### Sample Outputs

![title](sample.png)


These figures were generated using graph-tool (https://graph-tool.skewed.de/).
