This version of the IIRABM is described in: Cockrell, Chase, and Gary An. "Sepsis reconsidered: identifying novel metrics for behavioral landscape characterization with a high-performance computing implementation of an agent-based model." Journal of theoretical biology 430 (2017): 157-168; and Cockrell, Robert Chase, and Gary An. "Examining the controllability of sepsis using genetic algorithms on an agent-based model of systemic inflammation." PLoS computational biology 14, no. 2 (2018): e1005876.

Simulation parameters are all set in the “Parameters.h” file prior to compilation.  The IIRABM can be compiled using any C++ compiler that has the C++11 standard and access to an MPI library.  A typical compilation and execution would look like:

$CC -o IIRABM_executable *.cpp -std=c++11
$mpirun -n<Number of Processing Cores> IIRABM_executable


All modifiable parameters are defined in the “Parameters.h” file

Copyright (C) 2017  Gary An and Chase Cockrell

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
