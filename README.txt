---------------------------------------------------------------

                           Set up
			   
---------------------------------------------------------------

This program calculates the eigen-energies and functions using a
harmonic oscillator basis. NOTE, you need to have gsl installed on
your system.  Integrals for this program are computed using
gauss-hermite quadrature using a modified program which can be found
at http://www.mymathlib.com/quadrature/gauss_hermite.html.  I was not
able to find the author name for credit.


This program assumes you know:


  1.  You know the potential V(x)

  2.  the minimum x_m such that V'(x_m) = 0 and V(x_m +- dx) >= V(x_m)
  [rigt now it can only be zero]

The program calculates the eigenenergies and eigensates of a
system. The potential should only have one minima (I should expand it
at some point).


The optimum frequency is found by minimizing:

< 0 | H | 0 > 

where | 0 > is the lowest order harmonic oscillator eigenstate with
frquency w.  The frequency w which minimizes the expectation is the
frequency chosen.  If there is an error in obtaining w, the program
will exit and prompt the user manually enter w.

The parameters you set at the begining of QEQUICK.cc are

  1. me : mass of the particle.
  2. hbar: the value of the reduced planck constant.
  3. N : the number of basis functions (the dimension of the matrix)
  4. double V(double x): return the value of the potential at position "x".
  

To check convergence, calculate the eigen-energies for increasing
values of N.  When there is no change in the energies for increasing
N, you have convergence.

WARNING, WARNING: So far this program can only go as high as N=100 and
even lower for very fast growing potentials.  For N>100 the integrals
are prone to error and calculations of gamma functions become
intractible.


-----------------------------------------------------------------------

                         Running The Program

-----------------------------------------------------------------------

Once you've edited the program run one of the following commands:

   gcc -static QEQUICK.cc -o QEQUICK -lgsl -lgslcblas -lm
   gcc QEQUICK.cc -o QEQUICK -lgsl -lgslcblas -lm

depending on how you have set up gsl.  If it doesn't work, check your
gsl installation. Afterward run:

   ./QEQUICK

If you did not modify the main() program you should get a print-out of
the frequency and all the eigen-values.
