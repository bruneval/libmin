/*  This is my library */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


typedef struct { 
  int status;            // info
  int ndim;              // dimension of the minimization space
  int history_record;    // nhist
  double gtol;
  int iter;
  double *diag;
  double *work;

  int    line_info;
  int    line_infoc;
  int    line_nfev;
  double line_stp;
  double line_stpmin;
  double line_stpmax;
  double line_dginit;
  double line_finit;
  double line_stx;
  double line_fx;
  double line_dgx;
  double line_sty;
  double line_fy;
  double line_dgy;
  double line_stmin;
  double line_stmax;
  bool   line_bracket;
  bool   line_stage1;

} libmin_plan;


libmin_plan * libmin_init(int ndim, int history_record);
libmin_plan * libmin_init_diag(int ndim, int history_record, double *diag);

int libmin_execute(libmin_plan *p, double *x, double f, double *gradf);
void libmin_destroy(libmin_plan * p);
void lbfgs(int ndim, int history_record, double *x, double f, double *gradf, double *diag, double *work, int *status,
       double *gtol, double line_stpmin, double line_stpmax, double *line_stp, int *iter, int *line_info,
       int *line_nfev,
       double *line_dginit, double *line_finit,
       double *line_stx,  double *line_fx,  double *line_dgx,
       double *line_sty,  double *line_fy,  double *line_dgy,
       double *line_stmin,  double *line_stmax,
       bool *line_bracket, bool *line_stage1, int *line_infoc);


