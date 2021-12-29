SSDFO performs an approximate unconstrained global minimization of a not
necessarily smooth function of many continuous arguments. 
No gradients are needed.
A limited amount of noise is tolerated. 

The program was written by Morteza Kimiaei (University of Vienna). 
Please inform the author at kimiaeim83@univie.ac.at if you make 
serious use of this code. 

This directory contains 
* the Matlab source code for Version 3.0 (December 29, 2021), 
* a driver program driverSSDFO.m showing its use, 
* and the paper

M. Kimiaei, A. Neumaier, and P. Faramarzi
New subspace method for unconstrained black box optimization,
submitted (2021).


This paper describes the method implemented and some test results. 
Please cite this paper when using this package in scientific work.



SSDFO includes 

* FFD estimating the gradient by the forward finite difference
* FFDD estimating the directional derivative
* WolfeSearch  (getAccept, globalMin, goodBrack, goodStep, and
  interpolatCubic) finding a step size satisfying Wolfe condition,
* subspaceDir computing the subspace direction,
* updateSY updating the subspace information.

We rewrote and improved WolfeSearch according to lineSearch.m -- 
the optimization toolbox of Matlab software.

The old name was SSBBO. 
