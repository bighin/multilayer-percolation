# Monte Carlo simulation of percolation in multilayer systems

This code studies bond percolation in bilayer and generic multilayer systems by means of Monte Carlo simulation. The aim is studying the phase diagram a function of p -- the probability of activating an in-plane bond -- and p_perp -- the probability of activating a link between 

# Requirements

A modern C compiler (GCC or Clang) and the following libraries:

- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- NCurses

For compiling the code, the CMake build system is needed.

# Usage

The code uses the CMake system, therefore it will be compiled with the commands `mkdir build`, `cd build`, `cmake ..`.

Then go back to the main folder, and run the code with the command `./build/multilayer`.
