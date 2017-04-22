#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<malloc.h>

double pi = 3.14159265359;


double Gaussian(double x, double sigma );

// Function prototype for Runge-Kutta method for updated the wavefunction in time.
// RK3 is going to take the current wavefunction (double complex *wfn) at time
// t_i and update to its value at a future time t_f = t_i + dt 
// Arguments to the function:
//    arg1 - (int) dim - number of points in the wavefunction
//    arg2 - (double) *xvec - vector of x-values that the wfn is evaluated at
//    arg3 - (double complex) *wfn - vector that stores the wfn values at time t_i
//           when the function is called, and when the function has been executed,
//           this vector will hold the wfn values at time t_f = t_i + dt
//    arg4 - (double) dx - differential step along x between subsequent points in *xvec
//    arg5 - (double) dt - differential step along time between subsequent times (i.e. t_f - t_i)
void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);

// Function prototype for Second Derivative
// dfdt is going to take the wave function, take the second derivative 
// at each point multiplied by the imaginary unit 
// and return the time-derivative of the wave function. 
// arguments: arg 1 - dim, number of points in the wf
//            arg 2 - psivec - array of wavefunction values
//            arg 3 - dpsi - array of d/dt wavefunction values
//            arg 4  the difference between dx , increment along x axis 
void dfdt(int dim, double complex *psivec, double complex *dpsi, double dx );
double Collapse(int dim, double *xvec, double *P, double complex *wfn);

int main () {
       // evaluate wavefunction at this many points
       int dim = 1000;
       double *x, *P;
       double complex *wfn, *dpsi, *wfn_star;
       double *Cr, *Ci;
       // Length of domain in atomic units
       double L = 300.;
       // full width at half max of gaussin is ~2.35*sigma
       double sigma = 5;

	double dx;
	int i;

	dx = L/dim;

        x = (double *)malloc(dim*sizeof(double));
        wfn = (double complex *)malloc(dim*sizeof(double complex));
        dpsi= (double complex *)malloc(dim*sizeof(double complex));
        P   = (double *)malloc(dim*sizeof(double));
        wfn_star = (double complex *)malloc(dim*sizeof(double complex));
	
//Steve G Edits

	Cr = (double *)malloc(dim*sizeof(double));
	Ci = (double *)malloc(dim*sizeof(double));

//edits end

        for (i=0; i<dim; i++) {

          x[i] = (i-dim/2)*dx;

          wfn[i] = Gaussian(x[i], sigma) + 0.*I;
       //   printf("  %e  %e  %e\n",x[i],creal(wfn[i]),cimag(wfn[i]));
  }

  //RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);
  int pidx=1;
  for (int j=0; j<1; j++) {
  
        RK3(dim, x, wfn, dx, 0.001);
	//dfdt( dim, wfn, dpsi, dx);  
        if (j%50==0) {
          printf("\n\n#%i\n",pidx);
	  for (i=0; i<dim; i++) { 
            wfn_star[i] = conj(wfn[i]); // complex conjugate of wfn[i]
            P[i] = creal(wfn_star[i]*wfn[i]);
	   // printf(" %f %e %e %e\n",x[i],creal(wfn[i]),cimag(wfn[i]),creal(P[i]));

//Steve G Edits

            //Cr[i] = creal(wfn[i]);
            //Ci[i] = cimag(wfn[i]);
	    //printf(" %f %e %e %e\n",x[i],Cr[i],Ci[i]);
            //Collapse(dim, x, P, wfn);  
//edits end

          }
          Collapse(dim, x, P, wfn);
        pidx++;
        }   
}

//Collapse BS


}

void dfdt(int dim, double complex *psivec, double complex  *dpsi, double dx ) {
// write code here 
 int j;
 dpsi[0] = 0. + 0.*I;
 dpsi[dim] = 0. + 0.*I;

for (j=1; j<dim; j++) {
     
	dpsi[j] = (I/2.)*(psivec[j+1] - 2*psivec[j] + psivec[j-1])/(dx*dx);
}







}

void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt) {

  int i;
  double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;
  // Temporary arrays for computing derivatives of wfns and approximate updates to wfns
  wfn_dot = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn2 = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn3 = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn_np1 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k1 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k2 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k3 = (double complex *)malloc((dim+1)*sizeof(double complex));

  // Must initialize all (real and imaginary parts of) these elements of the arrays to zero 
  for (i=0; i<=dim; i++) {
    wfn_dot[i] = 0. + 0.*I;
    wfn2[i] = 0. + 0.*I;
    wfn3[i] = 0. + 0.*I;
    wfn_np1[i] = 0. + 0.*I;
    k1[i] = 0. + 0.*I;
    k2[i] = 0. + 0.*I;
    k3[i] = 0. + 0.*I;

  }

  // Get dPsi(n)/dt at initial time!	
  dfdt(dim, wfn, wfn_dot, dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k1[i] = dt*wfn_dot[i];
    wfn2[i] = wfn[i] + k1[i]/2.;
  }
  // Get dPsi(n+k1/2)/dt
  dfdt(dim, wfn2, wfn_dot, dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k2[i] = dt*wfn_dot[i];
    wfn3[i] = wfn[i] + k2[i]/2.;
  }
  // Get dPsi(n+k2/2)/dt
  dfdt(dim, wfn3, wfn_dot, dx);
  // Compute approximate update with Euler step
  // Then update actual wfn
  for (i=0; i<=dim; i++) {
    k3[i] = dt*wfn_dot[i];
    wfn_np1[i] = wfn[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
    wfn[i] = wfn_np1[i];
  }
  // wfn vector has now been updated!  
  // Now free memory associated with temporary vectors
  free(wfn_dot);
  free(wfn2);
  free(wfn3);
  free(wfn_np1);
  free(k1);
  free(k2);
  free(k3);

}



double Gaussian(double x, double sigma) {

  double prefac = 1./(sigma*sqrt(2*pi));
   double fun = exp(-x*x/(2*sigma*sigma));


  return prefac*fun;

}


double Collapse(int dim, double *xvec, double *P, double complex *wfn) {

  // TO DO:
  // (1) Find maximum value of P array

	int index;
	int length = dim; 
	double max = P[0];
        printf("  max is %12.10e\n",max);
        int *F;
        F = (int *)malloc(dim*sizeof(int));
	for(int i=1;i<length;i++){
		if(P[i]>max){
			max = P[i];
			index = i;
		}
	}
 
  // (2) Determine the frequency associated with this maximum value...
 
	int freq = 0;
  
	int max_freq = 1000;
	int freq_scale = floor(max_freq/max); 
        printf("  MAX IS %12.10e\n",max);
        printf("  FREQ SCALE %i\n", freq_scale);


	int Length2 = 0;
	for (int i=0; i<length ; i++){
	F[i] = floor(freq_scale*P[i]);
	Length2 = Length2 + F[i];
	printf(" freq %i \n", F[i]);
	

		
	}

        printf("  TOTAL LENGTH IS %i\n",Length2);

        



	//free(F); 
	






  // e.g. if max(P) -> 1e-3, maybe you want this to correspond to a frequency 
  // of 1000, in which case, you scale every other P(x) value by 1000.
  // (3) Construct a list of possible position values with each position repeated 
  //     a number of times N which is equal to floor( P(x)*1e6)... this will be called the choice array
  //     C[] and it will have M values.
	int *C;
	C = (int *)malloc(Length2*sizeof(int)); 
	int choice_count = 0;
	int maxj = 0;
	for (int i=0;i<length;i++){
		maxj = F[i];
		for (int j=0;j<maxj;j++){
			C[choice_count] = xvec[i];
			choice_count++;
		}
	}

	for (int d;d<Length2;d++)
		printf("  CHOICE ARRAY %i\n",C[d]);
 // (4) Choose a random integer i between 0 and M where M is the length of your choice array (C[]).
  //     The collapsed position will correspond to C[i]
	int range = Length2 - 1;
	int rchoice = rand() % range;
	double rand_xvec = C[rchoice];

  // (5) recalculate your wavefunction so that it corresponds to a Gaussian function centered
  //     at position C[i].

	free(F);
	free(C);
}


// END 
