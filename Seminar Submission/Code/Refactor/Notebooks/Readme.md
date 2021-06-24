---
title: HHL
tags: []
---

# Documentation to Advanced Quantum Computing Seminar Implementation

## Run

In order to run these notebooks, it is recommended to install Pluto Notebooks, as they have been tested on the same. 

    Pkg>add("Pluto"); 
    julia> import Pluto 
    julia> Pluto.run()


## HHL.jl 

Contains the main algortihm of the problem. For the sake of ease, the script generates a unit Poisson problem of arbitatry dimensions. You need only specify **nqubits**. 

Also note that the default number of qubits is set to neig = 12 for accuracy, atleast in the case of the Poisson Problem. 

Returns the HHL solution and the standard Julia Linear Algebra "\" solution. 

Main interface is driver. Usage is as follows:

    x_hhl,x_julia  = driver(nqubits,ndims=1)
    
I assume atleast for the sake of the submission that reading and plotting the solutions from this script is a trivial exercise. So an example is shown for the 1D case alone. 

## Controlled Rotation.jl 

Contains a vanilla bitwise to scalar implementation using successive Ry gates.

In addition it also contains a multiplexer circuit that is more versatile than the vanilla version of the implementation. 



## Phase Estimation.jl

Contains an implementation for the Phase Estimation circuit, referred from the Introduction to Quantum Computing Lecture at SCCS 5. 

Tested against the standard implementation and the eigenvalues are visualized and the values are then compared. 

## Quantum Fourier Transform.jl 

Naive implementation of QFT used inturn in Phase Estimation. Tested and benchmarked for some arbitary vectors and finally validated against a library implementation. 

## Example.jl 

Simple example that shows how the user interface be used.
