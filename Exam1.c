#include<stdio.h>
#include<stdlib.h> 
#include<math.h>

// calculating derivative pf psi_1(x)
int main() {


  int i, MAX_I;
  double x, xf, dx, fx, fpx, fxf;

  //rise = fxf-fx
  //run = xf - x
  //slipe = (fxf-fx)/(xf-x)

  double pi = 4*atan(1.0);
  //double pi = 3.141592653589793;
 
  MAX_I = 8;
  dx = 10./MAX_I;
  // psi_1(x) = sqrt(2/10)*sin(pi*x/10)
  for (i=0; i<=MAX_I; i++) { 

   // calculate numerical derivative 
   x = i*dx; 
   xf = (i+1)*dx;

   
   fx = sqrt(2/10.)*sin(pi*x/10.);
   fxf = sqrt(2/10.)*sin(pi*xf/10.);

  ror = (fxf-fx)/(xf-x); 
  
  //calculate analytical derivative 
  fpx = sqrt(2/10.)*(pi/10)*cos(pi*x/10);

  printf(" %i  %f  %f\n",i, x, ror, fpx);




 }




}
 


