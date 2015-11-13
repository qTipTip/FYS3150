# Project 1 - FYS3150
## Solving the Poisson Equation in one dimension

In the folder `/article/` is all the source code related to the report as well as the report itself.
In `/src/` is the source code used for computing the various numerical results discussed in the report. The programs store their data in the folder `/data/`.

Compile `main.cpp` with `g++ -larmadillo -std=c++11 main.cpp` and run
the framework with `python framework.py`. The python program will then perform all the calculations
and plotting for the various values of n.

The python framework does a bit too much at once, I should have sub-moduled it, but now it's one command for all the numerical solutions.

The graphical output is stored in `/plots/`.
