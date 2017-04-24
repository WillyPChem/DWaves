#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double pi = 4.0*atan(1.0);

// Function prototype
double Gaussian(double x, double sigma );

int main() {

  int i;
  double x[2000], g[2000];
  double gridspace = 0.1e-10;
  double sigma = 5e-9;

  for (i=0; i<2000; i++) {

    x[i] = (i-1000)*gridspace;
    // Evaluate Gaussian at x[i]
    g[i] = Gaussian(x[i], sigma);
    // Print xvalue and gaussian in scientific notation
    printf("  %e  %e\n",x[i],g[i]);
  }

}

// Function definition - defines what the function does
double Gaussian(double x, double sigma) {

  double prefac = 1./(sigma*sqrt(2*pi));
  double fun = exp(-x*x/(2*sigma*sigma));
 

  return prefac*fun;

}

