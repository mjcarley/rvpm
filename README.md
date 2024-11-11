RVPM is an implementation of the (Reformulated) Vortex Particle Method
based on the work of Alvarez and Ning:

https://doi.org/10.2514/1.J063045

The code is under development and liable to change, but it is possible
to do some useful things with it.

# Prerequisites

You need to have installed some other libraries, which are available
from the same place you got RVPM:
 - https://github.com/mjcarley/blaswrap
 - https://github.com/mjcarley/grbf
 - https://github.com/mjcarley/wbfmm
 - https://github.com/mjcarley/gqr

# Installation

RVPM uses GNU autotools for configuration and installation. The
following sequence should work on any reasonable system. Execute

`. autogen.sh`

in the root directory of the source tree. Then configure the code with

`./configure (options)`

and install with

`make`

and

`make install`

# Basic steps in solving a problem

Generate a set of particles, using `rvpm-grid` (or any other program
that suits you). For example,

`rvpm-grid -f vortex_ring -p 0.8 -p 0.1 -s 0.02 -a 1 > ring.dat`

uses a built-in function to approximate a vortex ring of radius 0.8
and core radius 0.1, using vortex particles of radius 0.02 with an
overlap of 1. 

Advance the solution using `rvpm-solve`:

`rvpm-solve -n 10 -d 0.05 < ring.dat  > solution.dat`

advances the solution for 10 time steps of size dt=0.05, placing the
final result in `solution.dat`.

Optionally, regrid the solution onto a new set of particles to restart
the time stepping:

`rvpm-grid -r solution.dat -a 1 -s 0.02 > restart.dat`
