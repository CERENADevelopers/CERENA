#ifndef _MY_geneExp_MCM3
#define _MY_geneExp_MCM3

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <udata.h>
#include <math.h>
#include <mex.h>

             void fu_geneExp_MCM3(void *user_data, double t);
             void fsu_geneExp_MCM3(void *user_data, double t);
             void fv_geneExp_MCM3(realtype t, N_Vector x, N_Vector dx, void *user_data);
             void dvdx_geneExp_MCM3(realtype t, N_Vector x, N_Vector dx, void *user_data);
             void dvddx_geneExp_MCM3(realtype t, N_Vector x, N_Vector dx, void *user_data);
             void dvdu_geneExp_MCM3(realtype t, N_Vector x, N_Vector dx, void *user_data);
             void dvdp_geneExp_MCM3(realtype t, N_Vector x, N_Vector dx, void *user_data, int ip);
             void fx0_geneExp_MCM3(N_Vector x0, void *user_data);
             void fdx0_geneExp_MCM3(N_Vector x0, N_Vector dx0, void *user_data);
             int DenseJacFn_geneExp_MCM3(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector F, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             int fsx_geneExp_MCM3(int nps, realtype t, N_Vector x, N_Vector dx, N_Vector F, N_Vector *sx, N_Vector *sdx, N_Vector *sF, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
             void fsx0_geneExp_MCM3(int ip, N_Vector sx0, void *user_data);
             void fsdx0_geneExp_MCM3(int ip, N_Vector sdx0, void *user_data);
             void fy_geneExp_MCM3(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x);
             void dfydp_geneExp_MCM3(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int nps, int nobs);
             void dfydx_geneExp_MCM3(double t, double *dydx, double *y, double *p, double *k, double *u, double *x);
             void fsy_geneExp_MCM3(double t, int nt, int it, int ip, int sp, int nx, int ny, double *sy, double *p, double *k, double *u, double *x, double *su, double *sx);


#endif /* _MY_geneExp_MCM3 */
