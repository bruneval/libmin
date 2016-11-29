/*  This is my library */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/*
 logical,parameter,private :: diagco=.FALSE.
 integer,parameter,private :: iprint(2) = [-1,0]
 !internal set-up of the LBFGS subroutine ...
 integer,protected  :: MP,LP
*/


typedef struct { 
  int status;            // info
  int ndim;              // dimension of the minimization space
  int history_record;    // nhist
  double tolerance;      // eps
  double gtol;
  double stpmin;
  double stpmax;
  int iter;
  double *diag;
  double *work;
  double stp;

  int    line_info;
  int    line_infoc;
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


libmin_plan * libmin_init(int ndim, int history_record, double tolerance);
libmin_plan * libmin_init_diag(int ndim, int history_record, double tolerance, double *diag);

int libmin_execute(libmin_plan *p, double *x, double f, double *gradf);
void libmin_destroy(libmin_plan * p);
void lbfgs(int ndim, int history_record, double *x, double f, double *gradf, bool diagco, double *diag, double tolerance, double *work, int *status,
       double *gtol, double stpmin, double stpmax, double *stp, int *iter, int *line_info,
       double *line_dginit, double *line_finit,
       double *line_stx,  double *line_fx,  double *line_dgx,
       double *line_sty,  double *line_fy,  double *line_dgy,
       double *line_stmin,  double *line_stmax,
       bool *line_bracket, bool *line_stage1, int *line_infoc);


