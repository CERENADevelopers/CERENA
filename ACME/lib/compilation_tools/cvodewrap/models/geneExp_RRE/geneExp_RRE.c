#include "geneExp_RRE.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <symbolic_functions.c>
#include <udata.h>
#include <math.h>
#include <mex.h>

#define pi 3.141592653589793


 int xdot_geneExp_RRE(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*4);
xdot_tmp[0] = p[2]*x_tmp[2] - 1.0*p[3]*x_tmp[0];
xdot_tmp[1] = p[1]*x_tmp[2] - 1.0*p[0]*x_tmp[1] - 1.0*k[0]*p[6]*x_tmp[1]*x_tmp[3];
xdot_tmp[2] = p[0]*x_tmp[1] - 1.0*p[1]*x_tmp[2] + k[0]*p[6]*x_tmp[1]*x_tmp[3];
xdot_tmp[3] = p[4]*x_tmp[0] - 1.0*p[5]*x_tmp[3];

  for (ix=0; ix<4; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_geneExp_RRE(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*20);

  for (ixB=0; ixB<20; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_geneExp_RRE(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
{
  int iyp;
  int ip;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  int np = *data->np;
  int ny = *data->ny;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *qBdot_tmp = N_VGetArrayPointer(qBdot);
  memset(qBdot_tmp,0,sizeof(double)*5*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }

  for (iyp=0; iyp<5*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_geneExp_RRE(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*4);
x0_tmp[0] = -(1.0*(p[9]*(k[3] - 1.0) - 1.0*k[3]*k[17]))/k[0];
x0_tmp[1] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[2] = (k[2]*k[16])/k[0];
x0_tmp[3] = (k[4]*k[18])/k[0];
  
  
  return;
}


 int Jv_geneExp_RRE(N_Vector v, N_Vector Jv, realtype t,
  	N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *v_tmp = N_VGetArrayPointer(v);
  double *Jv_tmp = N_VGetArrayPointer(Jv);
  memset(Jv_tmp,0,sizeof(double)*4);
Jv_tmp[0] = p[2]*v_tmp[2] - 1.0*p[3]*v_tmp[0];
Jv_tmp[1] = p[1]*v_tmp[2] - 1.0*v_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*v_tmp[3]*x_tmp[1];
Jv_tmp[2] = v_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*v_tmp[2] + k[0]*p[6]*v_tmp[3]*x_tmp[1];
Jv_tmp[3] = p[4]*v_tmp[0] - 1.0*p[5]*v_tmp[3];

  for (ix=0; ix<4; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_geneExp_RRE(N_Vector vB, N_Vector JvB, realtype t,
  	N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *vB_tmp = N_VGetArrayPointer(vB);
  double *JvB_tmp = N_VGetArrayPointer(JvB);
  memset(JvB_tmp,0,sizeof(double)*4);
JvB_tmp[0] = p[3]*vB_tmp[0] - 1.0*p[4]*vB_tmp[3];
JvB_tmp[1] = vB_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*vB_tmp[2]*(p[0] + k[0]*p[6]*x_tmp[3]);
JvB_tmp[2] = p[1]*vB_tmp[2] - 1.0*p[2]*vB_tmp[0] - 1.0*p[1]*vB_tmp[1];
JvB_tmp[3] = p[5]*vB_tmp[3] + k[0]*p[6]*vB_tmp[1]*x_tmp[1] - 1.0*k[0]*p[6]*vB_tmp[2]*x_tmp[1];

  for (ix=0; ix<4; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_geneExp_RRE(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_geneExp_RRE(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_geneExp_RRE(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
J->data[0] = -1.0*p[3];
J->data[3] = p[4];
J->data[5] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[3];
J->data[6] = p[0] + k[0]*p[6]*x_tmp[3];
J->data[8] = p[2];
J->data[9] = p[1];
J->data[10] = -1.0*p[1];
J->data[13] = -1.0*k[0]*p[6]*x_tmp[1];
J->data[14] = k[0]*p[6]*x_tmp[1];
J->data[15] = -1.0*p[5];

  for (iJ=0; iJ<16; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_geneExp_RRE(realtype t, N_Vector x,
  	N_Vector xdot, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(J);
  J->rowvals[0] = 0;
  J->rowvals[1] = 3;
  J->rowvals[2] = 1;
  J->rowvals[3] = 2;
  J->rowvals[4] = 0;
  J->rowvals[5] = 1;
  J->rowvals[6] = 2;
  J->rowvals[7] = 1;
  J->rowvals[8] = 2;
  J->rowvals[9] = 3;
  J->colptrs[0] = 0;
  J->colptrs[1] = 2;
  J->colptrs[2] = 4;
  J->colptrs[3] = 7;
  J->colptrs[4] = 10;
J->data[0] = -1.0*p[3];
J->data[1] = p[4];
J->data[2] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[3];
J->data[3] = p[0] + k[0]*p[6]*x_tmp[3];
J->data[4] = p[2];
J->data[5] = p[1];
J->data[6] = -1.0*p[1];
J->data[7] = -1.0*k[0]*p[6]*x_tmp[1];
J->data[8] = k[0]*p[6]*x_tmp[1];
J->data[9] = -1.0*p[5];
  return(0);
}


 int JBBand_geneExp_RRE(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_geneExp_RRE(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_geneExp_RRE(long int N, realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, DlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
JB->data[0] = p[3];
JB->data[2] = -1.0*p[2];
JB->data[5] = p[0] + k[0]*p[6]*x_tmp[3];
JB->data[6] = -1.0*p[1];
JB->data[7] = k[0]*p[6]*x_tmp[1];
JB->data[9] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[3];
JB->data[10] = p[1];
JB->data[11] = -1.0*k[0]*p[6]*x_tmp[1];
JB->data[12] = -1.0*p[4];
JB->data[15] = p[5];

  for (iJ=0; iJ<16; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_geneExp_RRE(realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, SlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(JB);
  JB->rowvals[0] = 0;
  JB->rowvals[1] = 2;
  JB->rowvals[2] = 1;
  JB->rowvals[3] = 2;
  JB->rowvals[4] = 3;
  JB->rowvals[5] = 1;
  JB->rowvals[6] = 2;
  JB->rowvals[7] = 3;
  JB->rowvals[8] = 0;
  JB->rowvals[9] = 3;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 2;
  JB->colptrs[2] = 5;
  JB->colptrs[3] = 8;
  JB->colptrs[4] = 10;
  return(0);
}


 int sx_geneExp_RRE(int Ns, realtype t, N_Vector x, N_Vector xdot,
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double *sxdot_tmp = N_VGetArrayPointer(sxdot);
  memset(sxdot_tmp,0,sizeof(double)*4);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + x_tmp[1] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 1: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) + x_tmp[2] - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] - 1.0*x_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 2: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0] + x_tmp[2];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 3: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0] - 1.0*x_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 4: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3] + x_tmp[0];

  } break;

  case 5: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3] - 1.0*x_tmp[3];

  } break;

  case 6: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*x_tmp[1]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*x_tmp[1]*x_tmp[3] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 7: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 8: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  case 9: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[1]*sx_tmp[2] - 1.0*sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[2] = sx_tmp[1]*(p[0] + k[0]*p[6]*x_tmp[3]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[3]*x_tmp[1];
sxdot_tmp[3] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[3];

  } break;

  }
 for (ix=0; ix<4; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_geneExp_RRE(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*4);
  switch (ip) {
  case 9: {
sx0_tmp[0] = -(1.0*(k[3] - 1.0))/k[0];

  } break;

  }

  return;
}


void y_geneExp_RRE(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*1];
y[it+nt*1] = x[it+nt*2];
y[it+nt*2] = x[it+nt*0];
y[it+nt*3] = x[it+nt*3];
y[it+nt*4] = p[7]*x[it+nt*3] + p[8];
    
    return;
}


void dydp_geneExp_RRE(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  memset(dydp,0,sizeof(double)*5*np*nt);
  switch (plist[ip]) {
  case 7: {
dydp[it+nt*(4+ip*ny)] = x[it+nt*3];

  } break;

  case 8: {
dydp[it+nt*(4+ip*ny)] = 1.0;

  } break;

  }
  }
  
  return;
}


void dydx_geneExp_RRE(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*20);
dydx[2] = 1.0;
dydx[5] = 1.0;
dydx[11] = 1.0;
dydx[18] = 1.0;
dydx[19] = p[7];
  
  return;
}


void sy_geneExp_RRE(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)] + x[it+nt*3];

  } break;

  case 8: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)] + 1.0;

  } break;

  case 9: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(4+np*ny)] = p[7]*sx[it+nt*(3+np*nx)];

  } break;

  }
  
  return;
}
int root_geneExp_RRE(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_geneExp_RRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dr_dp;
  switch (ip) {
  }
  return(dr_dp);
}
double s2root_geneExp_RRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddr_dpdp;
  switch (ip) {
  }
  return(ddr_dpdp);
}
double srootval_geneExp_RRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dg_dp;
  switch (ip) {
  }
  return(dg_dp);
}
double s2rootval_geneExp_RRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddg_dpdp;
  switch (ip) {
  }
  return(ddg_dpdp);
}
