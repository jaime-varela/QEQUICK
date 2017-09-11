#ifndef MATH_FUNC_H
#define MATH_FUNC_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sum.h>


/************ MATH FUNCTIONS ****************/

const double pi = 3.14159265358979323846;


// subtracts matrix arrays
//A<-A-B

void matrix_sub(complex double *A, complex double *B, int N){
  
  int i,j;

  for(i=0; i < N; i++){
    for(j=0; j < N; j++){
      A[i+N*j] = A[i+N*j]-B[i+N*j];
    }
  }
}


//returns (I)^(n) where n>0.
complex double ipwr(int n) {
  if(n==0)
    return 1;
  else if (n==1)
    return I;
  else{
  n = abs(n);
  complex double fk;
  fk = (n%2==0)? 1: I;
  if(n%2==0)
    fk *= ((n/2)%2==0)? -1:1;
  else
    fk *= (((n-1)/2)%2==0)? -1:1;

  return fk;}
}

//sign function
int sgn(double x) {
  if(x>0) 
    return 1;
  else 
    return -1;
}

double Sign(double x) {
  if(x<0) 
    return -1.0;
  else 
    return 1.0;
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

double fact(int n) 
{
  if(n < 2)
    return 1;
  return n * fact(n-1);
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


/* il(x) = Sqrt(pi / 2*x) I_{l+1/2}(x) */
void mod_besseli(int l, double x, double *il, double *il_deriv) {
  double expx = exp(x);
  *il = expx * gsl_sf_bessel_il_scaled(l, x);
  *il_deriv = expx * (l * gsl_sf_bessel_il_scaled(l-1, x) +
		      (l+1) * gsl_sf_bessel_il_scaled(l+1, x)) / (2*l+1);
  return;
}

/* kl(x) = Sqrt(2 / pi*x) K_{l+1/2}(x) */
void mod_besselk(int l, double x, double *kl, double *kl_deriv) {
  /* there is a bug in the gsl normalization */
  double norm = 2/pi; 
  double expx = exp(-x) * norm;
  *kl = expx * gsl_sf_bessel_kl_scaled(l, x);
  *kl_deriv = -expx * (l * gsl_sf_bessel_kl_scaled(l-1, x) +
		       (l+1) * gsl_sf_bessel_kl_scaled(l+1, x)) / (2*l+1);
  return;
}

void get_real_imag(gsl_matrix_complex *A, gsl_matrix *B, int N, int k){
  int i,j;
  gsl_complex z;
  double v;
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      z = gsl_matrix_complex_get(A, i, j);
      v = (k==0)? z.dat[0] : z.dat[1];
      gsl_matrix_set(B, i, j, v);
    }
  }
}


/* kl(x) array */
void mod_besselk_array(int l, double x, double *kl) {
 /* there is a bug in the gsl normalization */
  double norm = 2/pi;
  double expx = exp(-x) * norm;
  gsl_sf_bessel_kl_scaled_array(l, x, kl);
  int j;
  for(j=0; j <= l; j++) 
    kl[j] *= expx;
}


/* pm^l(x): as defined by the routines in
   http://ceta.mit.edu/ceta/comp_spec_func/ which allow |x| > 1
   analytic continuation. 
*/

/* This function only returns the absolute value of the legendre
   polynomial.  The true value is Pm^l(x)=(-I)^m pm^l(x). 

*/
void legendre(int l, int m, double x, double *plm, double *plm_deriv) 
{
  double norm;
  norm = (m < 0)? fact(l-fabs(m))/fact(l+fabs(m)) : 1; // comes from definition
  m = fabs(m); 
  /* The complex and sign factors must be taken care of externally
    with the factor (-I)^m */

  int i,j;
  int M = m+1, 
    N = l+2,
    NM = N*M;

  double *p = (double*) malloc(sizeof(double) * NM),
    *pd = (double*) malloc(sizeof(double) * NM);

  for(j=0; j < M; j++) {
    for(i=0; i < N; i++) {      
      int index = j*N+i;
      p[index] = 0.0;
      pd[index] = 0.0;
    }
  }
  p[0] = 1.0;

  if(fabs(x) == 1.0) {
    double xi = x;
    for(i=1; i < N; i++) {
      p[i] = xi;
      xi *= x;
      pd[i] = 0.5*i*(i+1)*xi;
    }
    if(m > 0)
      for(i=1; i < N; i++) 
	pd[N + i] = 301;
    if(m > 1) {
      xi = x;
      for(i=1; i < N; i++) {
	xi *= x;
	pd[2*N+i] = -0.25 * (i+2)*(i+1)*i*(i-1)*xi; 
      }
    }  
  } else {   
    int ls=1;

    if(fabs(x) > 1.0)
      ls = -1;    

    double xq = sqrt(ls * (1.0 - x*x));
    double xs = ls * (1.0-x*x);
    
    for(i=1; i < M; i++) 
      p[i*N+i] = -ls * (2*i-1) * xq * p[(i-1)*N+(i-1)];
    for(i=0; i < M; i++) 
      p[i*N+i+1] = (2*i+1)*x*p[i*N+i];
    
    for(i=0; i < M; i++) {
      for(j=i+2; j < N; j++) 
	p[i*N+j] = ((2*j-1) * x * p[i*N+j-1] - (i+j-1)*p[i*N+j-2]) / (j-i);
    }
    pd[0] = 0.0;
    for(j=1; j < N; j++)
      pd[j] = ls * j * (p[j-1] - x * p[j]) / xs;
    for(i=1; i <= m; i++) {
      for(j=i; j <= l; j++)
	pd[i*N+j] = ls * i * x * p[i*N+j] / xs + (j+i) * (j-i+1) * p[(i-1)*N+j] / xq;	
    }
  }

  *plm = norm * fabs(p[m*N+l]);
  *plm_deriv = norm * fabs(pd[m*N+l]);  

  free(p);
  free(pd);

  return;
}

double wigner_3jsym(int l1, int l2, int l3, 
		int m1, int m2, int m3) {
  //selection rules
  if(abs(m1) > l1 || abs(m2) > l2 || abs(m3) > l3)
    return 0.0;
  else if(m1 + m2 != -m3)
    return 0.0;
  else if(l3 < abs(l1-l2) || l3 > l1 + l2)
    return 0.0;
  else if(m1 ==0 && m2==0 && m3==0 && (l1+l2+l3)%2==1)
    return 0.0;
  else
  return gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);

}

int blas_matrix_mult(gsl_matrix_complex *A, gsl_matrix_complex *B, 
		     gsl_matrix_complex *C, 
		     gsl_complex alpha, gsl_complex beta) {
  // C = alpha * AxB + beta * C
  return gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
			alpha, A, B,
			beta, C);
}

int blas_matrix_add(gsl_complex alpha, 
		    gsl_vector_complex *A, gsl_vector_complex *B) {
  // B = alpha * A + B
  return gsl_blas_zaxpy(alpha, A, B);
}

void matrix_mult(gsl_matrix_complex *A, gsl_matrix_complex *B, 
		 gsl_matrix_complex *C, int N) {
  // C = -AxB + C
  int i, j, k;
  for(i=0; i < N; i++) {
    for(j=0; j < N; j++) {
      gsl_complex sum, aik, bkj, cij;
      sum.dat[0] = sum.dat[1] = 0;
      for(k=0; k < N; k++) {
	aik = gsl_matrix_complex_get(A, i, k);
	bkj = gsl_matrix_complex_get(B, k, j);
	sum.dat[0] += aik.dat[0]*bkj.dat[0] - aik.dat[1]*bkj.dat[1];
	sum.dat[1] += aik.dat[0]*bkj.dat[1] + aik.dat[1]*bkj.dat[0];
      }
      cij = gsl_matrix_complex_get(C, i, j);
      sum.dat[0] = -sum.dat[0] + cij.dat[0];
      sum.dat[1] = -sum.dat[1] + cij.dat[1];
      gsl_matrix_complex_set(C, i, j, sum);
    }
  }
}

void transpose(gsl_matrix_complex *M, int N) {
  int i,j;
  gsl_complex mij, mji;

  for(i=0; i < N; i++) {
    for(j=i; j < N; j++) {
      mij = gsl_matrix_complex_get(M, i, j);
      mji = gsl_matrix_complex_get(M, j, i);
      gsl_matrix_complex_set(M, i, j, mji);
      gsl_matrix_complex_set(M, j, i, mij);
    }
  }
}
#endif
