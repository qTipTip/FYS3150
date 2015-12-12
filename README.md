# FYS3150 - Computational Physics

## Introduction 

This repository is ment to contain various files related to the course
Computational Physics (FYS3150) held at the University of Oslo, fall 2015.

I will hopefully spend some time updating this readme as more work is done in
the course, and hopefully all code present in the repository should have clear
instructions on how to be ran.

## Course outline

The course has five mandatory projects where it is required to upload the
associated report as well as source code.

Every week there are a set of weekly exercises, that I will try to finish and
I'll upload the code writted or related material here.

### List of Projects

#### Project 1

We solve the Poisson equation in its one-dimensional form for a spherically symmetric system.
The problem reduces to that of a **tridiagonal matrix problem** and we implement the **tridiagonal matrix algorithm**
in order to solve this. We also discuss nuances of the numerical approximation used in order to discretize our problem.

#### Project 2

We examine the Schrödinger equation for two electrons in a harmonic oscillator well.
It can be shown that this can be formulated as an **eigenvalue problem**, and we solve this
problem using the **Jacobi method**. We consider the system for various choices of the different parameters
and analyze under what conditions our results are reasonable. 
**Bogus alert!** Me not having taken physics, I'm probably way off on some of my points.

#### Project 3

We consider three numerical integration methods for computing the **ground
state correlation energy** between two electrons in a hydrogen atom. Among the
methods considered are Guassian-Legendre quadrature, Gaussian-Laguerre
quadrature and Monte Carlo-integration.

#### Project 4

We look at the Ising model for magnetic systems with no external magnetic field. The metropolis algorithm is well suited for this kind of statistical analysis, and will be used in this project. We first derive some analytical results specific to the 2x2-case and compare these to our numerical computations. The code is parallellized using OpenMPI which lets you run the code on several clusters.

#### Project 5

We simulate argon using face centered cubic lattice, the Lennard-Jones potential and the Velocity Verlet time integrator.
The melting temperature of argon is measured to be somewhere around 300K. We also look at the energy conservation of Velocity Verlet and how it compares to the Euler-Cromer algorithm. 
