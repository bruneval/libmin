/* testing */


#include "libmin.h"

int main()
{
 libmin_plan * plan;
 int ndim=3;
 int nhist=5;
 int info;
 double tolerance;
 double *xx;
 double ff;
 double *gradf;
 int i;

 tolerance = 1.0e-7;

 xx    = malloc(ndim*sizeof(double));
 gradf = malloc(ndim*sizeof(double));

 plan = libmin_init(ndim,nhist,tolerance);

 ff = -1.0;
 for(i=0; i<ndim; i++) {
   xx[i]    = 0.01;
   gradf[i] = 0.02;
 }

 info = libmin_execute(plan,xx,ff,gradf);

 libmin_destroy(plan);

 free(xx);
 free(gradf);

 return 0;

}

