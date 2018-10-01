# OptimalOBB
Optimal Oriented Bounding Box algorithms

## Contents

This repository contains several Matlab scripts implementing algorithms for finding the optimal oriented bounding box of a set of points.

In `src/HYBBRID`, you will find the reference implementation for the HYBBRID algoritm described in https://dl.acm.org/citation.cfm?id=2019641

In `src/ORourke`, you will find an implementation of O'Rourke's algorithm describe in https://cs.smith.edu/~jorourke/Papers/MinVolBox.pdf

In `src/Barequet`, you will find an implementation of Barequet and Har-Peled's algorithm in doi:10.1006/jagm.2000.1127

Other methods are available in their own files, including PCA and random walk

Note that a Matlab implementation of Korsawe's algorithm is available on https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box
