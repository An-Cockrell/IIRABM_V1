Simulation parameters are all set in the “Parameters.h” file prior to compilation.  The IIRABM can be compiled using any C++ compiler that has the C++11 standard and access to an MPI library.  A typical compilation and execution would look like:

$CC -o IIRABM_executable *.cpp -std=c++11
$mpirun -n<Number of Processing Cores> IIRABM_executable


All modifiable parameters are defined in the “Parameters.h” file