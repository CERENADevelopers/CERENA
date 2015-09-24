#include "JakStat_RRE.h"
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


 int xdot_JakStat_RRE(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*9);
xdot_tmp[0] = -1.0*(x_tmp[0] - 1.0*x_tmp[1])*p[3];
xdot_tmp[1] = -1.0*(x_tmp[1] - 1.0*x_tmp[2])*p[3];
xdot_tmp[2] = -1.0*(x_tmp[2] - 1.0*x_tmp[3])*p[3];
xdot_tmp[3] = -1.0*(x_tmp[3] - 1.0*x_tmp[4])*p[3];
xdot_tmp[4] = -1.0*(x_tmp[4] - 2.0*x_tmp[5])*p[3];
xdot_tmp[5] = -(1.0*(p[3]*p[6]*x_tmp[5] - 1.0*p[2]*p[5]*x_tmp[8]))/p[6];
xdot_tmp[6] = (p[3]*p[6]*x_tmp[0] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[5]*p[17]*x_tmp[6])/p[5];
xdot_tmp[7] = -(1.0*(2.0*p[1]*p[5]*(pow(x_tmp[7],2)) - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[4]*p[17]*x_tmp[6]))/p[4];
xdot_tmp[8] = (p[1]*p[5]*(pow(x_tmp[7],2)) - 1.0*p[2]*p[4]*x_tmp[8])/p[4];

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_JakStat_RRE(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*99);

  for (ixB=0; ixB<99; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_JakStat_RRE(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
  memset(qBdot_tmp,0,sizeof(double)*11*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }

  for (iyp=0; iyp<11*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_JakStat_RRE(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*9);
x0_tmp[0] = (k[8]*k[62])/p[6];
x0_tmp[1] = (k[7]*k[61])/p[6];
x0_tmp[2] = (k[6]*k[60])/p[6];
x0_tmp[3] = (k[5]*k[59])/p[6];
x0_tmp[4] = (k[4]*k[58])/p[6];
x0_tmp[5] = (k[3]*k[57])/p[6];
x0_tmp[6] = (k[0]*k[54] - 1.0*p[4]*p[5]*(k[0] - 1.0))/p[5];
x0_tmp[7] = (k[1]*k[55])/p[5];
x0_tmp[8] = (k[2]*k[56])/p[5];
  
  
  return;
}


 int Jv_JakStat_RRE(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*9);
Jv_tmp[0] = p[3]*v_tmp[1] - 1.0*p[3]*v_tmp[0];
Jv_tmp[1] = p[3]*v_tmp[2] - 1.0*p[3]*v_tmp[1];
Jv_tmp[2] = p[3]*v_tmp[3] - 1.0*p[3]*v_tmp[2];
Jv_tmp[3] = p[3]*v_tmp[4] - 1.0*p[3]*v_tmp[3];
Jv_tmp[4] = 2.0*p[3]*v_tmp[5] - 1.0*p[3]*v_tmp[4];
Jv_tmp[5] = (p[2]*p[5]*v_tmp[8])/p[6] - 1.0*p[3]*v_tmp[5];
Jv_tmp[6] = (p[3]*p[6]*v_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*v_tmp[6];
Jv_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*v_tmp[6] - (4.0*p[1]*p[5]*v_tmp[7]*x_tmp[7])/p[4];
Jv_tmp[8] = (2.0*p[1]*p[5]*v_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*v_tmp[8];

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_JakStat_RRE(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*9);
JvB_tmp[0] = p[3]*vB_tmp[0] - (1.0*p[3]*p[6]*vB_tmp[6])/p[5];
JvB_tmp[1] = p[3]*vB_tmp[1] - 1.0*p[3]*vB_tmp[0];
JvB_tmp[2] = p[3]*vB_tmp[2] - 1.0*p[3]*vB_tmp[1];
JvB_tmp[3] = p[3]*vB_tmp[3] - 1.0*p[3]*vB_tmp[2];
JvB_tmp[4] = p[3]*vB_tmp[4] - 1.0*p[3]*vB_tmp[3];
JvB_tmp[5] = p[3]*vB_tmp[5] - 2.0*p[3]*vB_tmp[4];
JvB_tmp[6] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*vB_tmp[6] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*vB_tmp[7];
JvB_tmp[7] = (4.0*p[1]*p[5]*vB_tmp[7]*x_tmp[7])/p[4] - (2.0*p[1]*p[5]*vB_tmp[8]*x_tmp[7])/p[4];
JvB_tmp[8] = p[2]*vB_tmp[8] - (1.0*p[2]*p[5]*vB_tmp[5])/p[6];

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_JakStat_RRE(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_JakStat_RRE(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_JakStat_RRE(long int N, realtype t, N_Vector x,
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
J->data[6] = (p[3]*p[6])/p[5];
J->data[9] = p[3];
J->data[10] = -1.0*p[3];
J->data[19] = p[3];
J->data[20] = -1.0*p[3];
J->data[29] = p[3];
J->data[30] = -1.0*p[3];
J->data[39] = p[3];
J->data[40] = -1.0*p[3];
J->data[49] = 2.0*p[3];
J->data[50] = -1.0*p[3];
J->data[60] = -1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
J->data[61] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
J->data[70] = -(4.0*p[1]*p[5]*x_tmp[7])/p[4];
J->data[71] = (2.0*p[1]*p[5]*x_tmp[7])/p[4];
J->data[77] = (p[2]*p[5])/p[6];
J->data[80] = -1.0*p[2];

  for (iJ=0; iJ<81; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_JakStat_RRE(realtype t, N_Vector x,
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
  J->rowvals[1] = 6;
  J->rowvals[2] = 0;
  J->rowvals[3] = 1;
  J->rowvals[4] = 1;
  J->rowvals[5] = 2;
  J->rowvals[6] = 2;
  J->rowvals[7] = 3;
  J->rowvals[8] = 3;
  J->rowvals[9] = 4;
  J->rowvals[10] = 4;
  J->rowvals[11] = 5;
  J->rowvals[12] = 6;
  J->rowvals[13] = 7;
  J->rowvals[14] = 7;
  J->rowvals[15] = 8;
  J->rowvals[16] = 5;
  J->rowvals[17] = 8;
  J->colptrs[0] = 0;
  J->colptrs[1] = 2;
  J->colptrs[2] = 4;
  J->colptrs[3] = 6;
  J->colptrs[4] = 8;
  J->colptrs[5] = 10;
  J->colptrs[6] = 12;
  J->colptrs[7] = 14;
  J->colptrs[8] = 16;
  J->colptrs[9] = 18;
J->data[0] = -1.0*p[3];
J->data[1] = (p[3]*p[6])/p[5];
J->data[2] = p[3];
J->data[3] = -1.0*p[3];
J->data[4] = p[3];
J->data[5] = -1.0*p[3];
J->data[6] = p[3];
J->data[7] = -1.0*p[3];
J->data[8] = p[3];
J->data[9] = -1.0*p[3];
J->data[10] = 2.0*p[3];
J->data[11] = -1.0*p[3];
J->data[12] = -1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
J->data[13] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
J->data[14] = -(4.0*p[1]*p[5]*x_tmp[7])/p[4];
J->data[15] = (2.0*p[1]*p[5]*x_tmp[7])/p[4];
J->data[16] = (p[2]*p[5])/p[6];
J->data[17] = -1.0*p[2];
  return(0);
}


 int JBBand_JakStat_RRE(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_JakStat_RRE(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_JakStat_RRE(long int N, realtype t, N_Vector x,
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
JB->data[1] = -1.0*p[3];
JB->data[10] = p[3];
JB->data[11] = -1.0*p[3];
JB->data[20] = p[3];
JB->data[21] = -1.0*p[3];
JB->data[30] = p[3];
JB->data[31] = -1.0*p[3];
JB->data[40] = p[3];
JB->data[41] = -2.0*p[3];
JB->data[50] = p[3];
JB->data[53] = -(1.0*p[2]*p[5])/p[6];
JB->data[54] = -(1.0*p[3]*p[6])/p[5];
JB->data[60] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
JB->data[69] = -1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17];
JB->data[70] = (4.0*p[1]*p[5]*x_tmp[7])/p[4];
JB->data[79] = -(2.0*p[1]*p[5]*x_tmp[7])/p[4];
JB->data[80] = p[2];

  for (iJ=0; iJ<81; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_JakStat_RRE(realtype t, N_Vector x,
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
  JB->rowvals[1] = 1;
  JB->rowvals[2] = 1;
  JB->rowvals[3] = 2;
  JB->rowvals[4] = 2;
  JB->rowvals[5] = 3;
  JB->rowvals[6] = 3;
  JB->rowvals[7] = 4;
  JB->rowvals[8] = 4;
  JB->rowvals[9] = 5;
  JB->rowvals[10] = 5;
  JB->rowvals[11] = 8;
  JB->rowvals[12] = 0;
  JB->rowvals[13] = 6;
  JB->rowvals[14] = 6;
  JB->rowvals[15] = 7;
  JB->rowvals[16] = 7;
  JB->rowvals[17] = 8;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 2;
  JB->colptrs[2] = 4;
  JB->colptrs[3] = 6;
  JB->colptrs[4] = 8;
  JB->colptrs[5] = 10;
  JB->colptrs[6] = 12;
  JB->colptrs[7] = 14;
  JB->colptrs[8] = 16;
  JB->colptrs[9] = 18;
  return(0);
}


 int sx_JakStat_RRE(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*9);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[17]*x_tmp[6] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[17]*x_tmp[6] + (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 1: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (2.0*p[5]*(pow(x_tmp[7],2)))/p[4] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (p[5]*(pow(x_tmp[7],2)))/p[4] - 1.0*p[2]*sx_tmp[8] + (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];

  } break;

  case 2: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[5]*x_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5] + (p[2]*p[5]*sx_tmp[8])/p[6];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*x_tmp[8] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 3: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0] - 1.0*x_tmp[0] + x_tmp[1];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1] - 1.0*x_tmp[1] + x_tmp[2];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2] - 1.0*x_tmp[2] + x_tmp[3];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3] - 1.0*x_tmp[3] + x_tmp[4];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4] - 1.0*x_tmp[4] + 2.0*x_tmp[5];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*x_tmp[5] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[6]*x_tmp[0])/p[5] + (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 4: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (2.0*p[1]*p[5]*(pow(x_tmp[7],2)) - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[4]*p[17]*x_tmp[6])/(pow(p[4],2)) + (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] + ((pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*x_tmp[6])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - (1.0*(p[1]*p[5]*(pow(x_tmp[7],2)) - 1.0*p[2]*p[4]*x_tmp[8]))/(pow(p[4],2)) - (1.0*p[2]*x_tmp[8])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 5: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*x_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5] + (p[2]*p[5]*sx_tmp[8])/p[6];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - (1.0*(p[3]*p[6]*x_tmp[0] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[5]*p[17]*x_tmp[6]))/(pow(p[5],2)) - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*x_tmp[6])/p[5];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (2.0*p[1]*(pow(x_tmp[7],2)))/p[4] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (p[1]*(pow(x_tmp[7],2)))/p[4] - 1.0*p[2]*sx_tmp[8] + (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];

  } break;

  case 6: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[3]*p[6]*x_tmp[5] - 1.0*p[2]*p[5]*x_tmp[8])/(pow(p[6],2)) - 1.0*p[3]*sx_tmp[5] - (1.0*p[3]*x_tmp[5])/p[6] + (p[2]*p[5]*sx_tmp[8])/p[6];
sxdot_tmp[6] = (p[3]*x_tmp[0])/p[5] + (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 7: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 8: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 9: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 10: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 11: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 12: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 13: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 14: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 15: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 16: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] + (pow(t,3))*exp(-1.0*t*p[16])*p[0]*p[17]*x_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - 1.0*(pow(t,3))*exp(-1.0*t*p[16])*p[0]*p[17]*x_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  case 17: {
sxdot_tmp[0] = p[3]*sx_tmp[1] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[3]*sx_tmp[2] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = p[3]*sx_tmp[3] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = p[3]*sx_tmp[4] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[3]*sx_tmp[5] - 1.0*p[3]*sx_tmp[4];
sxdot_tmp[5] = (p[2]*p[5]*sx_tmp[8])/p[6] - 1.0*p[3]*sx_tmp[5];
sxdot_tmp[6] = (p[3]*p[6]*sx_tmp[0])/p[5] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*x_tmp[6] - 1.0*(pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6];
sxdot_tmp[7] = (pow(t,2))*exp(-1.0*t*p[16])*p[0]*x_tmp[6] + (pow(t,2))*exp(-1.0*t*p[16])*p[0]*p[17]*sx_tmp[6] - (4.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4];
sxdot_tmp[8] = (2.0*p[1]*p[5]*sx_tmp[7]*x_tmp[7])/p[4] - 1.0*p[2]*sx_tmp[8];

  } break;

  }
 for (ix=0; ix<9; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_JakStat_RRE(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*9);
  switch (ip) {
  case 4: {
sx0_tmp[6] = 1.0 - 1.0*k[0];

  } break;

  case 5: {
sx0_tmp[6] = - (1.0*(k[0]*k[54] - 1.0*p[4]*p[5]*(k[0] - 1.0)))/(pow(p[5],2)) - (1.0*p[4]*(k[0] - 1.0))/p[5];
sx0_tmp[7] = -(1.0*k[1]*k[55])/(pow(p[5],2));
sx0_tmp[8] = -(1.0*k[2]*k[56])/(pow(p[5],2));

  } break;

  case 6: {
sx0_tmp[0] = -(1.0*k[8]*k[62])/(pow(p[6],2));
sx0_tmp[1] = -(1.0*k[7]*k[61])/(pow(p[6],2));
sx0_tmp[2] = -(1.0*k[6]*k[60])/(pow(p[6],2));
sx0_tmp[3] = -(1.0*k[5]*k[59])/(pow(p[6],2));
sx0_tmp[4] = -(1.0*k[4]*k[58])/(pow(p[6],2));
sx0_tmp[5] = -(1.0*k[3]*k[57])/(pow(p[6],2));

  } break;

  }

  return;
}


void y_JakStat_RRE(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*6];
y[it+nt*1] = x[it+nt*7];
y[it+nt*2] = x[it+nt*8];
y[it+nt*3] = x[it+nt*5];
y[it+nt*4] = x[it+nt*4];
y[it+nt*5] = x[it+nt*3];
y[it+nt*6] = x[it+nt*2];
y[it+nt*7] = x[it+nt*1];
y[it+nt*8] = x[it+nt*0];
y[it+nt*9] = (p[9]*x[it+nt*7] + 2.0*p[9]*x[it+nt*8] + p[4]*p[7])/p[4];
y[it+nt*10] = (p[10]*x[it+nt*6] + p[10]*x[it+nt*7] + 2.0*p[10]*x[it+nt*8] + p[4]*p[8])/p[4];
    
    return;
}


void dydp_JakStat_RRE(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  memset(dydp,0,sizeof(double)*11*np*nt);
  switch (plist[ip]) {
  case 4: {
dydp[it+nt*(9+ip*ny)] = p[7]/p[4] - (1.0*(p[9]*x[it+nt*7] + 2.0*p[9]*x[it+nt*8] + p[4]*p[7]))/(pow(p[4],2));
dydp[it+nt*(10+ip*ny)] = p[8]/p[4] - (1.0*(p[10]*x[it+nt*6] + p[10]*x[it+nt*7] + 2.0*p[10]*x[it+nt*8] + p[4]*p[8]))/(pow(p[4],2));

  } break;

  case 7: {
dydp[it+nt*(9+ip*ny)] = 1.0;

  } break;

  case 8: {
dydp[it+nt*(10+ip*ny)] = 1.0;

  } break;

  case 9: {
dydp[it+nt*(9+ip*ny)] = (x[it+nt*7] + 2.0*x[it+nt*8])/p[4];

  } break;

  case 10: {
dydp[it+nt*(10+ip*ny)] = (x[it+nt*6] + x[it+nt*7] + 2.0*x[it+nt*8])/p[4];

  } break;

  }
  }
  
  return;
}


void dydx_JakStat_RRE(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*99);
dydx[8] = 1.0;
dydx[18] = 1.0;
dydx[28] = 1.0;
dydx[38] = 1.0;
dydx[48] = 1.0;
dydx[58] = 1.0;
dydx[66] = 1.0;
dydx[76] = p[10]/p[4];
dydx[78] = 1.0;
dydx[86] = p[9]/p[4];
dydx[87] = p[10]/p[4];
dydx[90] = 1.0;
dydx[97] = (2.0*p[9])/p[4];
dydx[98] = (2.0*p[10])/p[4];
  
  return;
}


void sy_JakStat_RRE(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = p[7]/p[4] - (1.0*(p[9]*x[it+nt*7] + 2.0*p[9]*x[it+nt*8] + p[4]*p[7]))/(pow(p[4],2)) + (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = p[8]/p[4] - (1.0*(p[10]*x[it+nt*6] + p[10]*x[it+nt*7] + 2.0*p[10]*x[it+nt*8] + p[4]*p[8]))/(pow(p[4],2)) + (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4] + 1.0;
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 8: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4] + 1.0;

  } break;

  case 9: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (x[it+nt*7] + 2.0*x[it+nt*8])/p[4] + (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 10: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (x[it+nt*6] + x[it+nt*7] + 2.0*x[it+nt*8])/p[4] + (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 11: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 12: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 13: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 14: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 15: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 16: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  case 17: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(9+np*ny)] = (p[9]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[9]*sx[it+nt*(8+np*nx)])/p[4];
sy[it+nt*(10+np*ny)] = (p[10]*sx[it+nt*(6+np*nx)])/p[4] + (p[10]*sx[it+nt*(7+np*nx)])/p[4] + (2.0*p[10]*sx[it+nt*(8+np*nx)])/p[4];

  } break;

  }
  
  return;
}
int root_JakStat_RRE(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_JakStat_RRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_JakStat_RRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_JakStat_RRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_JakStat_RRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
