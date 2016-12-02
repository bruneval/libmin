/* testing */


#include <stdio.h>
#include <math.h>
#include "libmin.h"

double eval_function(double *xx) {
  return  0.500 * pow( xx[0] - 1.0 , 2.0 )  
        + 0.125 * pow( xx[1] - 3.0 , 2.0 )  
        + 3.000 * pow( xx[2] + 2.0 , 2.0 ) ; 
}
void eval_gradient(double *xx, double *gradf) {
  gradf[0] = 1.0  * ( xx[0] - 1.0 );
  gradf[1] = 0.25 * ( xx[1] - 3.0 );
  gradf[2] = 6.0  * ( xx[2] + 2.0 );
}


int main()
{
 libmin_plan * plan;
 int const niter = 20;
 int const ndim = 3;
 int const nhist = 5;
 double const tolerance = 1.0e-7;
 double grad_norm2;
 double *xx;
 double *diag;
 double ff;
 double *gradf;
 int i;
 int iter;
 int info;


 xx    = malloc(ndim*sizeof(double));
 gradf = malloc(ndim*sizeof(double));
 diag  = malloc(ndim*sizeof(double));

 diag[0] = 1.0 / 1.00;
 diag[1] = 1.0 / 0.25;
 diag[2] = 1.0 / 6.00;

 plan = libmin_init_diag(ndim,nhist,diag);

 /* initial coordinates */
 xx[0] =  1.11258 + 1.5415 * 1;
 xx[1] =  1.11258 + 1.5415 * 2;
 xx[2] =  1.11258 + 1.5415 * 3;
 
 for(iter=0; iter<niter; iter++) {

   ff = eval_function(xx);
   eval_gradient(xx,gradf);

   printf(" --- iteration: %d  , f= %e \n",iter,ff);
   printf("X enter: %e %e %e \n",xx[0],xx[1],xx[2]);

   grad_norm2 = 0.0;
   for(i=0; i<ndim; i++) {
     grad_norm2 += gradf[i] * gradf[i];
   }

   if( grad_norm2 < tolerance ) {
     printf("Convergence reached\n");
     break;
   } 

   info = libmin_execute(plan,xx,ff,gradf);

   printf("X exit:  %e %e %e \n",xx[0],xx[1],xx[2]);

 }

 printf("Final X:  %e %e %e \n",xx[0],xx[1],xx[2]);

 libmin_destroy(plan);

 free(diag);
 free(xx);
 free(gradf);

 return 0;

}

