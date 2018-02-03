---------------------------------------------------------------

                           Introduction
			   
---------------------------------------------------------------

This program calculates the eigen-energies and functions using a
harmonic oscillator basis. So far it can calculate the 100x100 
Hamiltonian of a quartic potential in ~8e-3 seconds.

Integrals for this program are computed using gauss-hermite 
quadrature using a modified program which can be found
at http://www.mymathlib.com/quadrature/gauss_hermite.html.  I was not
able to find the author name for credit.

This program assumes you know:

  1.  You know the potential V(x)

  2.  the minimum x_m such that V'(x_m) = 0 and V(x_m +- dx) >= V(x_m)
  [rigt now it can only be zero, if it isn't just shift your coordinates.]

The program calculates the eigenenergies and eigensates of a
system. The potential should only have one minima (I should expand it
at some point).

The optimum frequency is found by minimizing:

< 0 | H | 0 > 

where | 0 > is the lowest order harmonic oscillator eigenstate with
frquency w.  The frequency w which minimizes the expectation is the
frequency chosen.  If there is an error in obtaining w, the program
will exit and prompt the user manually enter w.


---------------------------------------------------------------------------------------------

                                               Set up

---------------------------------------------------------------------------------------------

The parameters you set at the begining of QEQUICK.cpp are

  1. me : mass of the particle.
  2. hbar: the value of the reduced planck constant.
  3. iseven: a boolean that is true if the potential is even and odd otherwise.  

You define the potential of the system in the file "potential.h".  The potential takes in a double
and returns a double.  For example, if V(x) = x*x , replace the potential in "potential.h" with:

double V(double x){
    return x*x;
}

The definition of "V" is crucial and it can't be named anything else.

If you do not have GSL installed on your system then in the file, "math_functions.h", 
change the line:

   #define GSL_INSTALL 1
   
to be

   #define GSL_INSTALL 0


-----------------------------------------------------------------------

                         Running The Program

-----------------------------------------------------------------------

Once you've edited the program you can compile it with:

   g++ -std=c++11 -static QEQUICK.cpp -o QEQUICK -lgsl -lgslcblas -lm

if you have GSL, or:

   g++ -std=c++11 QEQUICK.cpp -o QEQUICK

if you do not have GSL. Depending on your system you may or may not need the
"-std=c++11" and "-static" flags.  

After compliation you can run:

   ./QEQUICK

If you don't have GSL installed the program will print the 40x40 approximate 
Hamiltonian as well as the frequency of the basis functions.  If you have GSL
installed the program will print the frequency and the eigenvalues of the 
approximate Hamiltonian.

You can change the amount of basis functions used in the approximation as well as
how many eigenvalues are printed in the command line.  For example:

   ./QEQUICK 60 20

Will use 60 basis functions and print 20 eigenvalues (if you don't have GSL then
the program will just print a 60x60 matrix).

You can also run the command:

  ./QEQUICK total
  
This command will generate the files "hamiltonian.dat","eigenvector.dat","eigenvalues.dat",
and "frequency.dat".  

The file "hamiltonian.dat" contains the numerical values of the 
hamiltonian in tab delimited form.  The jth line of "eigenvector.dat"
contains the jth eigenvector.  The value in the ith line and jth column of 
"eigenvalues.dat" contains the ith eigenenergy calculated 
using j+3 basis functions.  The file "frequency.dat" contains the frequency
of the basis functions used in the approximation.


WARNING: So far this program can only go as high as N=100 and
even lower for very fast growing potentials.  For N>100 the integrals
are prone to error and calculations of the hermite polynomials become
intractible.


-----------------------------------------------------------------------

                         Python programs

-----------------------------------------------------------------------

Also in this repo are some python programs that can be used to plot the
eigenenergies to check for convergence.
