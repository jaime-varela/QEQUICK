#ifndef MATH_FUNC_H
#define MATH_FUNC_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <vector>
#include <string>

#define GSL_INSTALL 1

#if GSL_INSTALL
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#endif

/************ MATH FUNCTIONS ****************/

const double pi = 3.14159265358979323846;
const double EXP= 2.71828182845904523536;
//sign function
int sgn(double x) {
  return (x>0)? 1:1;
}
double Sign(double x) {
  return (x<0)? -1.0:1.0;
}   
int max(int a, int b) {
  return (a < b)? b : a;
}
int kron(int a, int b) {
  return (a==b)? 1 : 0;
}

int isigma2(int a, int b){
  if (a==b)
    return 0;
  else if (b>a)
    return -1;
  else
    return 1;
}

//TODO: make this more efficient
double fact(int n) 
{
  if(n < 2)
    return 1;
  return n * fact(n-1);
}

bool isintodd(int a){
  if(a < 0)
    a = - a;
  return (a % 2 == 0)? false:true;
}

//standard bisection algorithim
double bisection(double error,double *params, double x0,double xf, 
double(*func)(double, void*))
{
  //  printf("%g \t %g \n",x0,xf);
  int sign = sgn(xf-x0);
  double yi = func(x0,params);
  double yf = func(xf,params);
  double ym = func(0.5*(xf+x0),params);
  double X0 = x0;
  double Xf = xf;
  double xm = 0.5*(xf+x0);
  //  printf("%g \n", xm);
  int i = 0;
  while(abs(Xf-X0)>error){
    if(sgn(yf)==sgn(ym)){
      Xf = xm;
      xm = 0.5*(Xf+X0);
      yf = ym;
      ym = func(xm,params);
    }
    else{
      X0 = xm;
      xm = 0.5*(Xf+X0);
      yi = ym;
      ym = func(xm,params);
    }
  }
  return xm;
}


//this function is usefull if you know there is a root somewhere
//AND it is close by
double findroot(double error,double *params, double x0,double step, 
double(*func)(double, void*))
{
  double y0=func(x0,params);
  int sign = (y0<0)? 1:-1;
  double y1=func(x0+sign*step,params);
  double y2=func(x0+sign*2*step,params);
  int i = 2*sign;
  if( y1>0 || y2>0){
    double XF = (y1 > 0)? (x0+step):(x0+2*step);
    return bisection(error,params,x0,XF,func);
  }
  else{
    while(sgn(y2) == -sign){
      i = i+ sign*1;
      y2 = func(x0+i*step,params);
    }
    //    printf("%g \n",x0+(i)*step);
    return bisection(error, params,x0+(i-1)*step,x0+i*step,func);
  }
}

//this function is a quadratic root finder
//use only if you know there is a root to the right
//of the initial point
double findroot3(double error,double *params, double x0,double step, 
double(*func)(double, void*))
{

  double y0=func(x0,params);
  int sign = (y0<0)? 1:-1;

  
  double x1 = x0+sign*step;
  double y1=func(x1,params);
  double x2 = x0+2*sign*step;
  double y2=func(x2,params);
  int i = 2*sign;
  int sign2 = sign;
  
  double a = 0;
  double b = 0;
  double c = 0;
  double root = 0;
  double disc = 0;
  while(sgn(y0) == sgn(y2)){
    a = pow((x0-x1)*(x0-x2)*(x1-x2),-1)*
      (x2*(y1-y0)+x1*(y0-y2)+x0*(y2-y1));
    b = pow((x0-x1)*(x0-x2)*(x1-x2),-1)*
      (x2*x2*(y0-y1)+x0*x0*(y1-y2)+x1*x1*(y2-y0));
    c = pow((x0-x1)*(x0-x2)*(x1-x2),-1)*
      (x2*(x1*(x1-x2)*y0+x0*(x2-x0)*y1)+x0*(x0-x1)*x1*y2);
    disc = b*b-4*a*c;
    
    if(disc < 0){
      double m = -1*(y2-y0);//note the negative sign
      double b = 0.5*(y2+y0)-m*x1;
      root = -b*pow(m,-1);
    }
    else{
      root = pow(2*a,-1)*(-b+sqrt(disc));
  }
    x0 = root;
    y0 = func(x0,params);
    sign2 = (y0 <0)? 1:-1;
    x1 = x0+sign*step;
    y1 = func(x1,params);
    x2 = x0+2*sign*step;
    y2 = func(x2,params);    }
  return bisection(error, params,x0,x2,func);
  
}

double sqrtlowfactorial(int n){
  double expon = 0.0;
  for(int i=1;i <= n;i++){
    expon += log(1.0*i);
  }
  return pow(EXP,expon/2.0);
}

//calculates sqrt(n!) 
// since we won't be going beyond n~200, we will not use an efficient
//factorial algorithm here 
double sqrtfact(int n){
#if GSL_INSTALL
    return sqrt(gsl_sf_gamma((double)n + 1.0));
#else
    return sqrtlowfactorial(n);
#endif
}

double Hermite_phys(int n, double x){
#if GSL_INSTALL
  return gsl_sf_hermite_phys(n,x);
#else
    double Hnm1 = 1.0;
    double Hn = 2.0 * x;
    if(n==0)
      return Hnm1;
    if(n==1)
      return Hn;
    double Hnp= Hn;
    double TWO = 2.0;
    double nn = 1.0*n;
    int k = 1;
    while(k<n){
      Hnp = TWO * (x * Hn - k * Hnm1);
      Hnm1 = Hn;
      Hn = Hnp;
      k++;
    }
    return Hnp;    
#endif
}


// -------------- Matrix routine (tested only for symmetric matrices) -------

//assumes the matrix is in an C-array
void print_matrix(double *A, int N){
  int i,j;
  printf("{");
  for(j=0;j<N;j++){
    printf("{");
    for(i=0;i<N;i++){
        if(i<N-1)
            printf("%g,",A[i+j*N]);
        else
            printf("%g",A[i+j*N]);
    }
    if(j < N-1){
        printf("},");
        printf("\n");
    }
    else{
        printf("}}\n");
    }
  }
}

std::string matrix_as_string(double* A, int N){
  std::ostringstream res;
  res.precision(17);
  int i,j;
  for(j=0;j<N;j++){
    for(i=0;i<N;i++){
        if(i<N-1){
            res << A[i+j*N];
            res << "\t";
        }
        else{
            res << A[i+j*N];            
        }
    }
    res << "\n";
    }
  return res.str();  
}

//populates Hsub_ij = H_ij  for 0 <= i,j < N
void populate_submatrix(double* H,double* Hsub,int N,int NTOT){
  int i,j;
  for(j=0;j<N;j++){
    for(i=0;i<N;i++){
      Hsub[i+j*N] = H[i+j*NTOT];
    }
  }
}

#if GSL_INSTALL
//the vectors E with the eigenvalues of the matrix H, Evec are the eigenvectors
// the column of Evec is the ith eigenvector.
void gsl_eigen_calculation(double* H,int N,double* E,double* Evec){

  //print_matrix(H,N);
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
    
    for (i = 0; i < N; i++) //right now only printing 40, change as desired.
      {
	double eval_i
	  = gsl_vector_get (eval, i);
	gsl_vector_view evec_i
	  = gsl_matrix_column(evec, i);	

	E[i] = eval_i;
  if(Evec != NULL){
    for(int ss=0; ss < N;ss++){
      //E[i] /= 1.0;
      //Evec[ss+]
      Evec[ss+i*N] = gsl_matrix_get(evec,ss,i);
    }
  }
      }
  }
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
}

std::string total_eigens_str(std::vector<std::vector<double> > & EE){
  std::ostringstream res;
  res.precision(17);
  int i,j;
  int N = EE.size();
//  std::cout << EE[0][12] << std::endl;
  for(j=0;j<N;j++){
    int nn = EE[j].size();
    for(i=0;i<nn;i++){
        if(i<nn-1){
            res << EE[j][i];
            res << "\t";
        }
        else{
            res << EE[j][i];            
        }
    }
    res << "\n";
    }
  return res.str();  
}
#endif

#endif
