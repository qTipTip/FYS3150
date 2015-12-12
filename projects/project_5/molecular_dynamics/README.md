# Atomic modelling of argon
## Project 5 - FYS3150

The project report can be found in the `article/` folder along with associated .tex files. 

In the `python_src/` folder you will find some auxilliary python scripts as well as the python framework used to run the simulations.
I made a class representing one single simulation of the system and used this to extract the relevant data from the simulations as needed.
Each self contained "task" was separated into its own python function. Rerunning the simulations is as simple as calling the relevant
function with a boolean flag (True/False) as argument. If true the system reruns the simulation, if not it uses old data. 

The C++ source code can be found in the `src/` folder, and was originally a fork of
https://github.com/andeplane/molecular-dynamics-fys3150
but some modifications has been made in accordance with the project problem text.

All relevant pictures can be found in the `pictures/` folder.
