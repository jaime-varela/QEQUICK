#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "math_functions.h"
#include "gauss_hermite_100pts_q.h"
#include "rootfinder.h"

#include "potential.h" //change this to whatever potential you have

//-------------------------------------------------------------
// Define your parameters here:
//-------------------------------------------------------------

//const double hbar = 1.054571800e-34;
//const double me = 9.10938e-31; // change the mass to whatever desired
const double hbar = 1.0;
const double me = 0.5;
const double coeff = hbar*hbar/(2*me);
const bool iseven = true;//set true if your potential is even, false otherwise
//const int N = 40;  //number of basis functions

/* The location of the minima (FOR NOW IT CAN ONLY BE ZERO). If your
 potential is not at zero, siply shift your coordinates eventually I
 will change it to allow for shifted potentials and multiple minima.*/
const double xm = 0.0;  

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
//const double W = 2.25116;

// --------------- Precompute V ---------------------------

const double xdiv = sqrt(me*W/hbar);
const int n100 = 100;
double *Vpvals = (double*)malloc((n100/2)*sizeof(double));
double *Vnvals = (double*)malloc((n100/2)*sizeof(double));
int populate_potential(double* Pos,double* Neg,int nn){
  for(int k = 0; k < nn/2;k++)
    Pos[k] = V((x[k])/xdiv);
  for(int k = 0; k < nn/2;k++)
    Neg[k] = V(-(x[k])/xdiv);
  return 0;
}
int vint = populate_potential(Vpvals,Vnvals,n100);
/*
  Here we define the harmonic oscillator basis functions WITHOUT the
  exponential factor.
*/
double Hbasis_func_no_exp(int v, double x){
  double c;
  double f;
  double one = 1.0;
  c = me*W/hbar;
  f = one/sqrt(pow(2,v));
  f *= one/sqrtfact(v);
  f *= pow(c/pi,0.25);
  double hermite = Hermite_phys(v, sqrt(c)*x);
  double result = hermite * f;
  return result;
}

//---------------- Precompute the basis functions --------------------------
double HDATAP[n100][n100/2];
double HDATAN[n100][n100/2];

int populate_hdata(double (&pos)[n100][n100/2],double (&neg)[n100][n100/2]){
  for(int i=0;i< n100;i++){
    for(int j=0;j<n100/2;j++){
      HDATAP[i][j] = Hbasis_func_no_exp(i,x[j]/xdiv);
      HDATAN[i][j] = Hbasis_func_no_exp(i,-x[j]/xdiv);//TODO: stop computing so many sqrts
    }
  }
  return 0;
}

int hdataint = populate_hdata(HDATAP,HDATAN);

/*
  The potential itegrand < i | x >  x) < x | j > with the appropriate variable
  transform x = u/c where c= sqrt(me W/hbar)
*/
double potential_element(double x, int i, int j){
  double c = sqrt(me*W/hbar);
  double result = (1/c)*Hbasis_func_no_exp(i, x/c)*V(x/c)*Hbasis_func_no_exp(j, x/c);
  return result;
}

double potential_element_mod(int xk, int i, int j,int sgn){
  double c = sqrt(me*W/hbar);
  double result = 0.0;
  if(sgn>0)
    result = (1/c)*HDATAP[i][xk]*Vpvals[xk]*HDATAP[j][xk];
  else
    result = (1/c)*HDATAN[i][xk]*Vnvals[xk]*HDATAN[j][xk];
  return result;
}

// the matrix elemens < i | H | j >
double matrix_element(int i, int j){
  double kinetic;
  if(i==j || i+2==j || i==j+2){
  kinetic = (hbar*W/4)*(kron(i,j)*sqrt(i*j) -
			       kron(i,j+2)*sqrt((j+1)*(j+2))+
			       kron(i+1,j+1)*sqrt((i+1)*(j+1))-
			       kron(i+2,j)*sqrt((i+1)*(i+2)));
  }
  else
    kinetic = 0.0;
  double result;
  //uses symmetry to reduce the number of integrands
  if(iseven && isintodd(i-j) )
    result = 0.0;
  else
    result = Gauss_Hermite_Integration_100pts_mod(i,j,potential_element_mod);
  result += kinetic;  
  return result;
}

//function which populates the matrix.
void populate_matrix(double *A, int N){
  int i,j;  //TODO: reduce calculation based on symmetry
  for(j=0;j<N;j++){
    for(i=0;i<N;i++){
      if(i<j)
        A[i+j*N] = A[j+i*N];
      else
        A[i+j*N] = matrix_element(i,j);
    }
  }
}


int main(int argc,char* argv[]){
  //if input is ./QEQUICK total generate data files. 
  if(argc==2 && !atoi(argv[1])){    
    std::string Sarg(argv[1]);
    if(Sarg != "total" && Sarg != "total " && Sarg != " total"){
      std::cout << "usage: \n./QEQUICK total" << std::endl;
      std::cout << "./QEQUICK [number of basis function] [number of eigenvals to print]" <<std::endl;
      std::cout << "example: ./QEQUICK 50 50" << std::endl;
    }
  std::ofstream hdata("hamiltonian.dat");  
  std::ofstream freq("frequency.dat");
  freq.precision(17);
  freq << W;
  freq.close();

  int ntot = 100;
  double* H  = (double*) malloc(ntot*ntot*sizeof(double)); //hamiltonian
  populate_matrix(H,ntot); //generating hamiltonian
  hdata << matrix_as_string(H,ntot);

#if GSL_INSTALL
  std::ofstream eigenval("eigenvalues.dat");
  std::ofstream eigenvec("eigenvector.dat");
  int nmin = 3;
  std::vector< std::vector<double> > EIGENVALUES(ntot);

  for(int ndata = nmin; ndata <= 100;ndata++){

    double* Hsub = (double*) malloc(ndata*ndata*sizeof(double));//subhamiltinion
    populate_submatrix(H,Hsub,ndata,ntot);
    //eigenvalue calculation
    double* Evals= (double*) malloc(ndata*sizeof(double));
    gsl_eigen_calculation(Hsub,ndata,Evals,NULL);
    
    for(int ss = 0; ss < ndata; ss++){
      EIGENVALUES[ss].push_back(Evals[ss]);
    }

    free(Evals);
    free(Hsub);
  }
  //std::cout << "here" <<std::endl;
  eigenval << total_eigens_str(EIGENVALUES);
  eigenval.close();

  double* Evec= (double*)malloc(ntot*ntot*sizeof(double));
  double* Evals= (double*)malloc(ntot*sizeof(double));
  gsl_eigen_calculation(H,ntot,Evals,Evec);

  eigenvec << matrix_as_string(Evec,ntot);
  eigenvec.close();
  free(Evec);
  free(Evals);
  free(H);
#else

  free(H);

#endif  
  
  hdata.close();

  }
  else{
  int N = 40; //default number of basis functions
  int Npr=40; // number of eigenvalues printed
  if(argc==2){
    N = atoi(argv[1]);
    Npr=N;
  }
  if(argc==3){
    N = atoi(argv[1]);
    Npr=atoi(argv[2]);
    if(Npr > N)
      Npr = N;
  }

  if(N > 100){
    printf("Basis must be smaller than 100 \n Setting N = 100 \n");
    N = 100;
    if(Npr > N)
      Npr = N;
  }
  double *H = (double*) malloc(N*N*sizeof(double)); //hamiltonian

  populate_matrix(H,N); //generating hamiltonian

#if GSL_INSTALL

  double *Evals = (double*) malloc(N*sizeof(double));
  gsl_eigen_calculation(H,N,Evals,NULL);
  int i;    
  for (i = 0; i < Npr; i++) //right now only printing 40, change as desired.
  	printf("eigenvalue %i = %g\n", i, Evals[i]);

  free(Evals);	
#else
  print_matrix(H,N);
#endif  

  //dealocations
  free(H);
  free(Vnvals);
  free(Vpvals);  
  }
  return 0;
}