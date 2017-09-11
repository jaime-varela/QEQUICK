#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "math_functions.h"
#include "gauss_hermite_100pts_q.h"
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "rootfinder.h"

//-------------------------------------------------------------
// Define your parameters here:
//-------------------------------------------------------------
//const double hbar = 1.054571800e-34;
//const double me = 9.10938e-31; // change the mass to whatever desired
const double hbar = 1.0;
const double me = 0.5;
const double coeff = hbar*hbar/(2*me);
const int N = 40;  //number of basis functions

/* The location of the minima (FOR NOW IT CAN ONLY BE ZERO). If your
 potential is not at zero, siply shift your coordinates eventually I
 will change it to allow for shifted potentials along with multiple
 minima.*/
const double xm = 0.0;  


//----------------------------------------------------------------
//----------------------- POTENTIAL ------------------------------
//----------------------------------------------------------------
/*
  Here you define the potential of the system V(x)

 */
double V(double x){
  double lambda = 0.1;
  return x*x + lambda*x*x*x*x;
}






//---------------------------------------------------------------
//--------------------- Defining W ------------------------------
//---------------------------------------------------------------




double expintegralofV(double u,double w){
  return (1-2*u*u)*V(u/sqrt(me*w/hbar)); 
}
    

double g(double w){
  return (2*sqrt(pi)*hbar*w/4) +
    Gauss_Hermite_Integration_100ptsW(w,expintegralofV);
}

double find_w(){
  double w0 = -(1/sqrt(2*exp(1)))*(V(1/(2*sqrt(2)) )
				  -V(sqrt(3/2.0) )); //initial guess
  double w;
  double step = 0.1;
  double err = 1e-3;

  w = rootfindbisection(g,w0,step,err);
  if(w < 0){ // TODO: isnan check
    printf("automatic method not viable, please supply w manually \n");
    exit (EXIT_FAILURE);
  }    
  else
    printf("W is given by : %g \n",w);
    return w;
  }

const double W = find_w();
//const double W = 10; // a user defined W



/*
  Here we define the harmonic oscillator basis functions WITHOUT the
  exponential factor.
*/
double Hbasis_func_no_exp(int v, double x){
  double c,f;
  c = me*W/hbar;
  f = 1/sqrt(pow(2,v)*gsl_sf_gamma((double)v + 1.0)); // TODO: fix the v_max = 170
  return f*pow(c/pi,0.25)*gsl_sf_hermite_phys(v, sqrt(c)*x);  
}





/*
  The potential itegrand < i | x >  V(x) < x | j > with the appropriate variable
  transform x = u/c where c= sqrt(me W/hbar)
*/
double potential_element(double x, int i, int j){
  double c = sqrt(me*W/hbar);
  double result = (1/c)*Hbasis_func_no_exp(i, x/c)*V(x/c)*Hbasis_func_no_exp(j, x/c);
  return result;
}


// the matrix elemens < i | H | j >
double matrix_element(int i, int j){
  
  double kinetic = (hbar*W/4)*(kron(i,j)*sqrt(i*j) -
			       kron(i,j+2)*sqrt((j+1)*(j+2))+
			       kron(i+1,j+1)*sqrt((i+1)*(j+1))-
			       kron(i+2,j)*sqrt((i+1)*(i+2)));
  double result;
  result = Gauss_Hermite_Integration_100pts(i,j,potential_element);
  result += kinetic;  
  return result;
}


//function which populates the matrix.
void populate_matrix(double *A, int N){
  int i,j;  //TODO: reduce calculation based on symmetry
  for(j=0;j<N;j++){
    for(i=0;i<N;i++){
      A[i+j*N] = matrix_element(i,j);
    }
  }
}




int main(){

  

  double *H = (double*) malloc(N*N*sizeof(double)); //hamiltonian
  

  populate_matrix(H,N); //generating hamiltonian


  /*
    The following calculates the eigenvalues and eigenstates for H.
    This program uses gsl routines to calculate eigenvalues.  If you
    prefer a different method just change the routines below.
   */
  gsl_matrix_view m
    = gsl_matrix_view_array (H, N, N);
  
  gsl_vector *eval = gsl_vector_alloc (N);
  gsl_matrix *evec = gsl_matrix_alloc (N, N);
  
  gsl_eigen_symmv_workspace * w =
    gsl_eigen_symmv_alloc (N);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  
  gsl_eigen_symmv_free (w);
  
  gsl_eigen_symmv_sort (eval, evec,
			GSL_EIGEN_SORT_ABS_ASC);
  
  {
    int i;
    
    for (i = 0; i < N; i++) //right now only printing 4, change as desired.
      {
	double eval_i
	  = gsl_vector_get (eval, i);
	gsl_vector_view evec_i
	  = gsl_matrix_column (evec, i);
	
	printf("eigenvalue %i = %g\n", i, eval_i);
	// remove "//" below if you want to print the eigenvectors.
	//printf ("eigenvector = \n");
	//gsl_vector_fprintf (stdout,
	//		    &evec_i.vector, "%g");
      }
  }
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  free(H);
  
  return 0;
}
