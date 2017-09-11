#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>






double max(double a, double b) {
  return ( a-b > 0)? a : b;
}


double min(double a, double b) {
  return ( a-b > 0)? b : a;
}


// returns the root given an inital guess x0.
// using the newton method.
//inputs are f a function, fd the derivative, x0 initial guess
// err is the error.  The truncation stops at |x_n+1 -x_n| < err. 
double rootfind(double (*f)(double)
		,double (*fd)(double)
		,double x0
		,double err) {

  double x1 = x0 - (f(x0)/fd(x0));
  while(fabs(x1-x0)> err){
    x0 = x1;
    x1 = x0 - (f(x0)/fd(x0));
  }


  return x1;
}



// uses the bisection method to find a root
double rootfindbisection(double (*f)(double)
			 ,double x0
			 ,double step
			 ,double err) {


  double y0 = f(x0);
  double y1 = f(x0 +step);

  double m = (y1-y0)/step;

  double x1 = -(((x0+step)*y0-x0*y1)/step)*(1/m);


  y1 = f(x1);
  int n = 3;
  while(y0/y1 > 0){
    y0 = y1;
    x0 = x1;
    y1 = f(x0 +step);
    m = (y1-y0)/step;
    x1 = -(((x0+step)*y0-x0*y1)/step)*(1/m);
    y1 = f(x1);
    n += 2;
  }


  double midpt;
  double left = min(x0,x1);
  double right = max(x0,x1);
  
  while(fabs(left-right)> err){
    midpt = (left+right)/2;
    if(f(midpt)/f(left) > 0){
      left = midpt;
    }
    else{
      right = midpt;
    }
    n += 2;
    
  }

  x1 = (left+right)/2;
  //printf("%d \n",n);
  return x1;
}



