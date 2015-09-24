#include "geneExp_EMRE.h"
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


 int xdot_geneExp_EMRE(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*18);
xdot_tmp[0] = p[2]*x_tmp[2] - 1.0*p[3]*x_tmp[0];
xdot_tmp[1] = p[4]*x_tmp[0] - 1.0*p[5]*x_tmp[1];
xdot_tmp[2] = p[0]*x_tmp[3] - 1.0*p[1]*x_tmp[2] + k[0]*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[1]*x_tmp[10] + k[0]*p[6]*x_tmp[3]*x_tmp[17];
xdot_tmp[3] = p[1]*x_tmp[2] - 1.0*p[0]*x_tmp[3] - 1.0*k[0]*p[6]*x_tmp[15] - 1.0*k[0]*p[6]*x_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*x_tmp[3]*x_tmp[17];
xdot_tmp[4] = (p[0]*x_tmp[10] + p[1]*x_tmp[11] - 2.0*k[0]*p[0]*x_tmp[4] + 2.0*k[0]*p[1]*x_tmp[12] - 2.0*(pow(k[0],2))*p[6]*x_tmp[4]*x_tmp[17] - 2.0*(pow(k[0],2))*p[6]*x_tmp[10]*x_tmp[15] + k[0]*p[6]*x_tmp[10]*x_tmp[17])/k[0];
xdot_tmp[5] = (p[0]*x_tmp[10] + p[1]*x_tmp[11] - 2.0*k[0]*p[1]*x_tmp[5] + 2.0*k[0]*p[0]*x_tmp[12] + 2.0*(pow(k[0],2))*p[6]*x_tmp[10]*x_tmp[16] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12]*x_tmp[17] + k[0]*p[6]*x_tmp[10]*x_tmp[17])/k[0];
xdot_tmp[6] = (p[3]*x_tmp[7] + p[2]*x_tmp[11] - 2.0*k[0]*p[3]*x_tmp[6] + 2.0*k[0]*p[2]*x_tmp[14])/k[0];
xdot_tmp[7] = p[2]*x_tmp[11] - 1.0*p[3]*x_tmp[7];
xdot_tmp[8] = (p[4]*x_tmp[7] + p[5]*x_tmp[17] + 2.0*k[0]*p[4]*x_tmp[9] - 2.0*k[0]*p[5]*x_tmp[8])/k[0];
xdot_tmp[9] = p[4]*x_tmp[6] - 1.0*p[3]*x_tmp[9] - 1.0*p[5]*x_tmp[9] + p[2]*x_tmp[16];
xdot_tmp[10] = p[1]*x_tmp[11] - 1.0*p[0]*x_tmp[10] - 1.0*k[0]*p[6]*x_tmp[10]*x_tmp[17];
xdot_tmp[11] = p[0]*x_tmp[10] - 1.0*p[1]*x_tmp[11] + k[0]*p[6]*x_tmp[10]*x_tmp[17];
xdot_tmp[12] = -(1.0*((pow(k[0],2))*p[1]*x_tmp[12] - 1.0*(pow(k[0],2))*p[1]*x_tmp[5] - 1.0*(pow(k[0],2))*x_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) + (pow(k[0],2))*x_tmp[12]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[0]*x_tmp[10] + k[0]*p[1]*x_tmp[11] - 1.0*(pow(k[0],3))*p[6]*x_tmp[10]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[10]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[10]*x_tmp[16]))/(pow(k[0],2));
xdot_tmp[13] = p[2]*x_tmp[12] - 1.0*p[0]*x_tmp[13] + p[1]*x_tmp[14] - 1.0*p[3]*x_tmp[13] - 1.0*k[0]*p[6]*x_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*x_tmp[13]*x_tmp[17];
xdot_tmp[14] = p[2]*x_tmp[5] + p[0]*x_tmp[13] - 1.0*p[1]*x_tmp[14] - 1.0*p[3]*x_tmp[14] + k[0]*p[6]*x_tmp[9]*x_tmp[10] + k[0]*p[6]*x_tmp[13]*x_tmp[17];
xdot_tmp[15] = p[1]*x_tmp[16] - 1.0*p[0]*x_tmp[15] + p[4]*x_tmp[13] - 1.0*p[5]*x_tmp[15] - 1.0*k[0]*p[6]*x_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*x_tmp[15]*x_tmp[17];
xdot_tmp[16] = p[0]*x_tmp[15] - 1.0*p[1]*x_tmp[16] + p[4]*x_tmp[14] - 1.0*p[5]*x_tmp[16] + k[0]*p[6]*x_tmp[8]*x_tmp[10] + k[0]*p[6]*x_tmp[15]*x_tmp[17];
xdot_tmp[17] = p[4]*x_tmp[7] - 1.0*p[5]*x_tmp[17];

  for (ix=0; ix<18; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_geneExp_EMRE(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*288);

  for (ixB=0; ixB<288; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_geneExp_EMRE(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
  memset(qBdot_tmp,0,sizeof(double)*16*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }

  for (iyp=0; iyp<16*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_geneExp_EMRE(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*18);
x0_tmp[7] = -(1.0*(p[9]*(k[3] - 1.0) - 1.0*k[3]*k[17]))/k[0];
x0_tmp[10] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[11] = (k[2]*k[16])/k[0];
x0_tmp[17] = (k[4]*k[18])/k[0];
  
  
  return;
}


 int Jv_geneExp_EMRE(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*18);
Jv_tmp[0] = p[2]*v_tmp[2] - 1.0*p[3]*v_tmp[0];
Jv_tmp[1] = p[4]*v_tmp[0] - 1.0*p[5]*v_tmp[1];
Jv_tmp[2] = v_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*v_tmp[2] + k[0]*p[6]*v_tmp[15] + k[0]*p[6]*v_tmp[1]*x_tmp[10] + k[0]*p[6]*v_tmp[10]*x_tmp[1] + k[0]*p[6]*v_tmp[17]*x_tmp[3];
Jv_tmp[3] = p[1]*v_tmp[2] - 1.0*v_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*v_tmp[15] - 1.0*k[0]*p[6]*v_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*v_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*v_tmp[17]*x_tmp[3];
Jv_tmp[4] = 2.0*p[1]*v_tmp[12] + (p[1]*v_tmp[11])/k[0] + (v_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*v_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*v_tmp[17])/k[0] - 2.0*k[0]*p[6]*v_tmp[15]*x_tmp[10];
Jv_tmp[5] = (p[1]*v_tmp[11])/k[0] - 2.0*p[1]*v_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*v_tmp[12])/k[0] + (v_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*v_tmp[17])/k[0] + 2.0*k[0]*p[6]*v_tmp[16]*x_tmp[10];
Jv_tmp[6] = 2.0*p[2]*v_tmp[14] - 2.0*p[3]*v_tmp[6] + (p[3]*v_tmp[7])/k[0] + (p[2]*v_tmp[11])/k[0];
Jv_tmp[7] = p[2]*v_tmp[11] - 1.0*p[3]*v_tmp[7];
Jv_tmp[8] = 2.0*p[4]*v_tmp[9] - 2.0*p[5]*v_tmp[8] + (p[4]*v_tmp[7])/k[0] + (p[5]*v_tmp[17])/k[0];
Jv_tmp[9] = p[4]*v_tmp[6] + p[2]*v_tmp[16] - 1.0*v_tmp[9]*(p[3] + p[5]);
Jv_tmp[10] = p[1]*v_tmp[11] - 1.0*v_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*v_tmp[17]*x_tmp[10];
Jv_tmp[11] = v_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*v_tmp[11] + k[0]*p[6]*v_tmp[17]*x_tmp[10];
Jv_tmp[12] = p[1]*v_tmp[5] + v_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*v_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*v_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*v_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*v_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*v_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*v_tmp[16]*x_tmp[10];
Jv_tmp[13] = p[2]*v_tmp[12] - 1.0*v_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*v_tmp[14] - 1.0*k[0]*p[6]*v_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*v_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*v_tmp[17]*x_tmp[13];
Jv_tmp[14] = p[2]*v_tmp[5] - 1.0*v_tmp[14]*(p[1] + p[3]) + v_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*v_tmp[9]*x_tmp[10] + k[0]*p[6]*v_tmp[10]*x_tmp[9] + k[0]*p[6]*v_tmp[17]*x_tmp[13];
Jv_tmp[15] = p[1]*v_tmp[16] - 1.0*v_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*v_tmp[13] - 1.0*k[0]*p[6]*v_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*v_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*v_tmp[17]*x_tmp[15];
Jv_tmp[16] = p[4]*v_tmp[14] - 1.0*v_tmp[16]*(p[1] + p[5]) + v_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*v_tmp[8]*x_tmp[10] + k[0]*p[6]*v_tmp[10]*x_tmp[8] + k[0]*p[6]*v_tmp[17]*x_tmp[15];
Jv_tmp[17] = p[4]*v_tmp[7] - 1.0*p[5]*v_tmp[17];

  for (ix=0; ix<18; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_geneExp_EMRE(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*18);
JvB_tmp[0] = p[3]*vB_tmp[0] - 1.0*p[4]*vB_tmp[1];
JvB_tmp[1] = p[5]*vB_tmp[1] - 1.0*k[0]*p[6]*vB_tmp[2]*x_tmp[10] + k[0]*p[6]*vB_tmp[3]*x_tmp[10];
JvB_tmp[2] = p[1]*vB_tmp[2] - 1.0*p[2]*vB_tmp[0] - 1.0*p[1]*vB_tmp[3];
JvB_tmp[3] = vB_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*vB_tmp[2]*(p[0] + k[0]*p[6]*x_tmp[17]);
JvB_tmp[4] = ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*vB_tmp[4])/k[0] - 1.0*vB_tmp[12]*(p[0] + k[0]*p[6]*x_tmp[17]);
JvB_tmp[5] = 2.0*p[1]*vB_tmp[5] - 1.0*p[1]*vB_tmp[12] - 1.0*p[2]*vB_tmp[14];
JvB_tmp[6] = 2.0*p[3]*vB_tmp[6] - 1.0*p[4]*vB_tmp[9];
JvB_tmp[7] = p[3]*vB_tmp[7] - 1.0*p[4]*vB_tmp[17] - (1.0*p[3]*vB_tmp[6])/k[0] - (1.0*p[4]*vB_tmp[8])/k[0];
JvB_tmp[8] = 2.0*p[5]*vB_tmp[8] + k[0]*p[6]*vB_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*vB_tmp[16]*x_tmp[10];
JvB_tmp[9] = vB_tmp[9]*(p[3] + p[5]) - 2.0*p[4]*vB_tmp[8] + k[0]*p[6]*vB_tmp[13]*x_tmp[10] - 1.0*k[0]*p[6]*vB_tmp[14]*x_tmp[10];
JvB_tmp[10] = vB_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*vB_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[17]) + (vB_tmp[12]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*vB_tmp[4]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*vB_tmp[5]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] - 1.0*k[0]*p[6]*vB_tmp[2]*x_tmp[1] + k[0]*p[6]*vB_tmp[3]*x_tmp[1] + k[0]*p[6]*vB_tmp[13]*x_tmp[9] - 1.0*k[0]*p[6]*vB_tmp[14]*x_tmp[9] + k[0]*p[6]*vB_tmp[15]*x_tmp[8] - 1.0*k[0]*p[6]*vB_tmp[16]*x_tmp[8];
JvB_tmp[11] = p[1]*vB_tmp[11] - 1.0*p[1]*vB_tmp[10] - 1.0*p[2]*vB_tmp[7] - (1.0*p[1]*vB_tmp[4])/k[0] - (1.0*p[1]*vB_tmp[5])/k[0] - (1.0*p[2]*vB_tmp[6])/k[0] + (p[1]*vB_tmp[12])/k[0];
JvB_tmp[12] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*vB_tmp[12])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[13] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*vB_tmp[5])/k[0] - 2.0*p[1]*vB_tmp[4];
JvB_tmp[13] = vB_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) - 1.0*p[4]*vB_tmp[15] - 1.0*vB_tmp[14]*(p[0] + k[0]*p[6]*x_tmp[17]);
JvB_tmp[14] = vB_tmp[14]*(p[1] + p[3]) - 1.0*p[1]*vB_tmp[13] - 1.0*p[4]*vB_tmp[16] - 2.0*p[2]*vB_tmp[6];
JvB_tmp[15] = vB_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) - 1.0*vB_tmp[16]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*vB_tmp[2] + k[0]*p[6]*vB_tmp[3] + 2.0*k[0]*p[6]*vB_tmp[4]*x_tmp[10] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[10];
JvB_tmp[16] = vB_tmp[16]*(p[1] + p[5]) - 1.0*p[1]*vB_tmp[15] - 1.0*p[2]*vB_tmp[9] - 2.0*k[0]*p[6]*vB_tmp[5]*x_tmp[10] + k[0]*p[6]*vB_tmp[12]*x_tmp[10];
JvB_tmp[17] = p[5]*vB_tmp[17] + (vB_tmp[12]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[5]*vB_tmp[8])/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*vB_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*vB_tmp[5])/k[0] - 1.0*k[0]*p[6]*vB_tmp[2]*x_tmp[3] + k[0]*p[6]*vB_tmp[3]*x_tmp[3] + k[0]*p[6]*vB_tmp[10]*x_tmp[10] - 1.0*k[0]*p[6]*vB_tmp[11]*x_tmp[10] + k[0]*p[6]*vB_tmp[13]*x_tmp[13] - 1.0*k[0]*p[6]*vB_tmp[14]*x_tmp[13] + k[0]*p[6]*vB_tmp[15]*x_tmp[15] - 1.0*k[0]*p[6]*vB_tmp[16]*x_tmp[15];

  for (ix=0; ix<18; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_geneExp_EMRE(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_geneExp_EMRE(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_geneExp_EMRE(long int N, realtype t, N_Vector x,
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
J->data[1] = p[4];
J->data[19] = -1.0*p[5];
J->data[20] = k[0]*p[6]*x_tmp[10];
J->data[21] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[36] = p[2];
J->data[38] = -1.0*p[1];
J->data[39] = p[1];
J->data[56] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[57] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[76] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17]))/k[0];
J->data[84] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[95] = -2.0*p[1];
J->data[102] = p[1];
J->data[104] = p[2];
J->data[114] = -2.0*p[3];
J->data[117] = p[4];
J->data[132] = p[3]/k[0];
J->data[133] = -1.0*p[3];
J->data[134] = p[4]/k[0];
J->data[143] = p[4];
J->data[152] = -2.0*p[5];
J->data[159] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[160] = k[0]*p[6]*x_tmp[10];
J->data[170] = 2.0*p[4];
J->data[171] = - 1.0*p[3] - 1.0*p[5];
J->data[175] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[176] = k[0]*p[6]*x_tmp[10];
J->data[182] = k[0]*p[6]*x_tmp[1];
J->data[183] = -1.0*k[0]*p[6]*x_tmp[1];
J->data[184] = (p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17])/k[0];
J->data[185] = (p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17])/k[0];
J->data[190] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[191] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[192] = -(1.0*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2));
J->data[193] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[194] = k[0]*p[6]*x_tmp[9];
J->data[195] = -1.0*k[0]*p[6]*x_tmp[8];
J->data[196] = k[0]*p[6]*x_tmp[8];
J->data[202] = p[1]/k[0];
J->data[203] = p[1]/k[0];
J->data[204] = p[2]/k[0];
J->data[205] = p[2];
J->data[208] = p[1];
J->data[209] = -1.0*p[1];
J->data[210] = -(1.0*p[1])/k[0];
J->data[220] = 2.0*p[1];
J->data[221] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])/k[0];
J->data[228] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17])))/(pow(k[0],2));
J->data[229] = p[2];
J->data[247] = - 1.0*p[0] - 1.0*p[3] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[248] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[249] = p[4];
J->data[258] = 2.0*p[2];
J->data[265] = p[1];
J->data[266] = - 1.0*p[1] - 1.0*p[3];
J->data[268] = p[4];
J->data[272] = k[0]*p[6];
J->data[273] = -1.0*k[0]*p[6];
J->data[274] = -2.0*k[0]*p[6]*x_tmp[10];
J->data[282] = k[0]*p[6]*x_tmp[10];
J->data[285] = - 1.0*p[0] - 1.0*p[5] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[286] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[293] = 2.0*k[0]*p[6]*x_tmp[10];
J->data[297] = p[2];
J->data[300] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[303] = p[1];
J->data[304] = - 1.0*p[1] - 1.0*p[5];
J->data[308] = k[0]*p[6]*x_tmp[3];
J->data[309] = -1.0*k[0]*p[6]*x_tmp[3];
J->data[310] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10]))/k[0];
J->data[311] = (2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])/k[0];
J->data[314] = p[5]/k[0];
J->data[316] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[317] = k[0]*p[6]*x_tmp[10];
J->data[318] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2));
J->data[319] = -1.0*k[0]*p[6]*x_tmp[13];
J->data[320] = k[0]*p[6]*x_tmp[13];
J->data[321] = -1.0*k[0]*p[6]*x_tmp[15];
J->data[322] = k[0]*p[6]*x_tmp[15];
J->data[323] = -1.0*p[5];

  for (iJ=0; iJ<324; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_geneExp_EMRE(realtype t, N_Vector x,
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
  J->rowvals[1] = 1;
  J->rowvals[2] = 1;
  J->rowvals[3] = 2;
  J->rowvals[4] = 3;
  J->rowvals[5] = 0;
  J->rowvals[6] = 2;
  J->rowvals[7] = 3;
  J->rowvals[8] = 2;
  J->rowvals[9] = 3;
  J->rowvals[10] = 4;
  J->rowvals[11] = 12;
  J->rowvals[12] = 5;
  J->rowvals[13] = 12;
  J->rowvals[14] = 14;
  J->rowvals[15] = 6;
  J->rowvals[16] = 9;
  J->rowvals[17] = 6;
  J->rowvals[18] = 7;
  J->rowvals[19] = 8;
  J->rowvals[20] = 17;
  J->rowvals[21] = 8;
  J->rowvals[22] = 15;
  J->rowvals[23] = 16;
  J->rowvals[24] = 8;
  J->rowvals[25] = 9;
  J->rowvals[26] = 13;
  J->rowvals[27] = 14;
  J->rowvals[28] = 2;
  J->rowvals[29] = 3;
  J->rowvals[30] = 4;
  J->rowvals[31] = 5;
  J->rowvals[32] = 10;
  J->rowvals[33] = 11;
  J->rowvals[34] = 12;
  J->rowvals[35] = 13;
  J->rowvals[36] = 14;
  J->rowvals[37] = 15;
  J->rowvals[38] = 16;
  J->rowvals[39] = 4;
  J->rowvals[40] = 5;
  J->rowvals[41] = 6;
  J->rowvals[42] = 7;
  J->rowvals[43] = 10;
  J->rowvals[44] = 11;
  J->rowvals[45] = 12;
  J->rowvals[46] = 4;
  J->rowvals[47] = 5;
  J->rowvals[48] = 12;
  J->rowvals[49] = 13;
  J->rowvals[50] = 13;
  J->rowvals[51] = 14;
  J->rowvals[52] = 15;
  J->rowvals[53] = 6;
  J->rowvals[54] = 13;
  J->rowvals[55] = 14;
  J->rowvals[56] = 16;
  J->rowvals[57] = 2;
  J->rowvals[58] = 3;
  J->rowvals[59] = 4;
  J->rowvals[60] = 12;
  J->rowvals[61] = 15;
  J->rowvals[62] = 16;
  J->rowvals[63] = 5;
  J->rowvals[64] = 9;
  J->rowvals[65] = 12;
  J->rowvals[66] = 15;
  J->rowvals[67] = 16;
  J->rowvals[68] = 2;
  J->rowvals[69] = 3;
  J->rowvals[70] = 4;
  J->rowvals[71] = 5;
  J->rowvals[72] = 8;
  J->rowvals[73] = 10;
  J->rowvals[74] = 11;
  J->rowvals[75] = 12;
  J->rowvals[76] = 13;
  J->rowvals[77] = 14;
  J->rowvals[78] = 15;
  J->rowvals[79] = 16;
  J->rowvals[80] = 17;
  J->colptrs[0] = 0;
  J->colptrs[1] = 2;
  J->colptrs[2] = 5;
  J->colptrs[3] = 8;
  J->colptrs[4] = 10;
  J->colptrs[5] = 12;
  J->colptrs[6] = 15;
  J->colptrs[7] = 17;
  J->colptrs[8] = 21;
  J->colptrs[9] = 24;
  J->colptrs[10] = 28;
  J->colptrs[11] = 39;
  J->colptrs[12] = 46;
  J->colptrs[13] = 50;
  J->colptrs[14] = 53;
  J->colptrs[15] = 57;
  J->colptrs[16] = 63;
  J->colptrs[17] = 68;
  J->colptrs[18] = 81;
J->data[0] = -1.0*p[3];
J->data[1] = p[4];
J->data[2] = -1.0*p[5];
J->data[3] = k[0]*p[6]*x_tmp[10];
J->data[4] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[5] = p[2];
J->data[6] = -1.0*p[1];
J->data[7] = p[1];
J->data[8] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[9] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[10] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17]))/k[0];
J->data[11] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[12] = -2.0*p[1];
J->data[13] = p[1];
J->data[14] = p[2];
J->data[15] = -2.0*p[3];
J->data[16] = p[4];
J->data[17] = p[3]/k[0];
J->data[18] = -1.0*p[3];
J->data[19] = p[4]/k[0];
J->data[20] = p[4];
J->data[21] = -2.0*p[5];
J->data[22] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[23] = k[0]*p[6]*x_tmp[10];
J->data[24] = 2.0*p[4];
J->data[25] = - 1.0*p[3] - 1.0*p[5];
J->data[26] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[27] = k[0]*p[6]*x_tmp[10];
J->data[28] = k[0]*p[6]*x_tmp[1];
J->data[29] = -1.0*k[0]*p[6]*x_tmp[1];
J->data[30] = (p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17])/k[0];
J->data[31] = (p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17])/k[0];
J->data[32] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[33] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[34] = -(1.0*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2));
J->data[35] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[36] = k[0]*p[6]*x_tmp[9];
J->data[37] = -1.0*k[0]*p[6]*x_tmp[8];
J->data[38] = k[0]*p[6]*x_tmp[8];
J->data[39] = p[1]/k[0];
J->data[40] = p[1]/k[0];
J->data[41] = p[2]/k[0];
J->data[42] = p[2];
J->data[43] = p[1];
J->data[44] = -1.0*p[1];
J->data[45] = -(1.0*p[1])/k[0];
J->data[46] = 2.0*p[1];
J->data[47] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])/k[0];
J->data[48] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17])))/(pow(k[0],2));
J->data[49] = p[2];
J->data[50] = - 1.0*p[0] - 1.0*p[3] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[51] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[52] = p[4];
J->data[53] = 2.0*p[2];
J->data[54] = p[1];
J->data[55] = - 1.0*p[1] - 1.0*p[3];
J->data[56] = p[4];
J->data[57] = k[0]*p[6];
J->data[58] = -1.0*k[0]*p[6];
J->data[59] = -2.0*k[0]*p[6]*x_tmp[10];
J->data[60] = k[0]*p[6]*x_tmp[10];
J->data[61] = - 1.0*p[0] - 1.0*p[5] - 1.0*k[0]*p[6]*x_tmp[17];
J->data[62] = p[0] + k[0]*p[6]*x_tmp[17];
J->data[63] = 2.0*k[0]*p[6]*x_tmp[10];
J->data[64] = p[2];
J->data[65] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[66] = p[1];
J->data[67] = - 1.0*p[1] - 1.0*p[5];
J->data[68] = k[0]*p[6]*x_tmp[3];
J->data[69] = -1.0*k[0]*p[6]*x_tmp[3];
J->data[70] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10]))/k[0];
J->data[71] = (2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])/k[0];
J->data[72] = p[5]/k[0];
J->data[73] = -1.0*k[0]*p[6]*x_tmp[10];
J->data[74] = k[0]*p[6]*x_tmp[10];
J->data[75] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2));
J->data[76] = -1.0*k[0]*p[6]*x_tmp[13];
J->data[77] = k[0]*p[6]*x_tmp[13];
J->data[78] = -1.0*k[0]*p[6]*x_tmp[15];
J->data[79] = k[0]*p[6]*x_tmp[15];
J->data[80] = -1.0*p[5];
  return(0);
}


 int JBBand_geneExp_EMRE(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_geneExp_EMRE(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_geneExp_EMRE(long int N, realtype t, N_Vector x,
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
JB->data[18] = -1.0*p[4];
JB->data[19] = p[5];
JB->data[37] = -1.0*k[0]*p[6]*x_tmp[10];
JB->data[38] = p[1];
JB->data[39] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
JB->data[46] = -1.0*k[0]*p[6]*x_tmp[1];
JB->data[51] = -1.0*k[0]*p[6];
JB->data[53] = -1.0*k[0]*p[6]*x_tmp[3];
JB->data[55] = k[0]*p[6]*x_tmp[10];
JB->data[56] = -1.0*p[1];
JB->data[57] = p[0] + k[0]*p[6]*x_tmp[17];
JB->data[64] = k[0]*p[6]*x_tmp[1];
JB->data[69] = k[0]*p[6];
JB->data[71] = k[0]*p[6]*x_tmp[3];
JB->data[76] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])/k[0];
JB->data[82] = -(1.0*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0];
JB->data[83] = -(1.0*p[1])/k[0];
JB->data[84] = -2.0*p[1];
JB->data[87] = 2.0*k[0]*p[6]*x_tmp[10];
JB->data[89] = (2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])/k[0];
JB->data[95] = 2.0*p[1];
JB->data[100] = -(1.0*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0];
JB->data[101] = -(1.0*p[1])/k[0];
JB->data[102] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17]))/k[0];
JB->data[106] = -2.0*k[0]*p[6]*x_tmp[10];
JB->data[107] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10]))/k[0];
JB->data[114] = 2.0*p[3];
JB->data[115] = -(1.0*p[3])/k[0];
JB->data[119] = -(1.0*p[2])/k[0];
JB->data[122] = -2.0*p[2];
JB->data[133] = p[3];
JB->data[137] = -1.0*p[2];
JB->data[151] = -(1.0*p[4])/k[0];
JB->data[152] = 2.0*p[5];
JB->data[153] = -2.0*p[4];
JB->data[161] = -(1.0*p[5])/k[0];
JB->data[168] = -1.0*p[4];
JB->data[171] = p[3] + p[5];
JB->data[178] = -1.0*p[2];
JB->data[190] = p[0] + k[0]*p[6]*x_tmp[17];
JB->data[191] = -1.0*p[1];
JB->data[197] = k[0]*p[6]*x_tmp[10];
JB->data[208] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
JB->data[209] = p[1];
JB->data[215] = -1.0*k[0]*p[6]*x_tmp[10];
JB->data[220] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
JB->data[221] = -1.0*p[1];
JB->data[226] = (k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16])/(pow(k[0],2));
JB->data[227] = p[1]/k[0];
JB->data[228] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))/(pow(k[0],2));
JB->data[231] = -1.0*k[0]*p[6]*x_tmp[10];
JB->data[232] = k[0]*p[6]*x_tmp[10];
JB->data[233] = ((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12])/(pow(k[0],2));
JB->data[243] = k[0]*p[6]*x_tmp[10];
JB->data[244] = k[0]*p[6]*x_tmp[9];
JB->data[246] = -1.0*p[2];
JB->data[247] = p[0] + p[3] + k[0]*p[6]*x_tmp[17];
JB->data[248] = -1.0*p[1];
JB->data[251] = k[0]*p[6]*x_tmp[13];
JB->data[257] = -1.0*p[2];
JB->data[261] = -1.0*k[0]*p[6]*x_tmp[10];
JB->data[262] = -1.0*k[0]*p[6]*x_tmp[9];
JB->data[265] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
JB->data[266] = p[1] + p[3];
JB->data[269] = -1.0*k[0]*p[6]*x_tmp[13];
JB->data[278] = k[0]*p[6]*x_tmp[10];
JB->data[280] = k[0]*p[6]*x_tmp[8];
JB->data[283] = -1.0*p[4];
JB->data[285] = p[0] + p[5] + k[0]*p[6]*x_tmp[17];
JB->data[286] = -1.0*p[1];
JB->data[287] = k[0]*p[6]*x_tmp[15];
JB->data[296] = -1.0*k[0]*p[6]*x_tmp[10];
JB->data[298] = -1.0*k[0]*p[6]*x_tmp[8];
JB->data[302] = -1.0*p[4];
JB->data[303] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[17];
JB->data[304] = p[1] + p[5];
JB->data[305] = -1.0*k[0]*p[6]*x_tmp[15];
JB->data[313] = -1.0*p[4];
JB->data[323] = p[5];

  for (iJ=0; iJ<324; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_geneExp_EMRE(realtype t, N_Vector x,
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
  JB->rowvals[2] = 0;
  JB->rowvals[3] = 1;
  JB->rowvals[4] = 1;
  JB->rowvals[5] = 2;
  JB->rowvals[6] = 3;
  JB->rowvals[7] = 10;
  JB->rowvals[8] = 15;
  JB->rowvals[9] = 17;
  JB->rowvals[10] = 1;
  JB->rowvals[11] = 2;
  JB->rowvals[12] = 3;
  JB->rowvals[13] = 10;
  JB->rowvals[14] = 15;
  JB->rowvals[15] = 17;
  JB->rowvals[16] = 4;
  JB->rowvals[17] = 10;
  JB->rowvals[18] = 11;
  JB->rowvals[19] = 12;
  JB->rowvals[20] = 15;
  JB->rowvals[21] = 17;
  JB->rowvals[22] = 5;
  JB->rowvals[23] = 10;
  JB->rowvals[24] = 11;
  JB->rowvals[25] = 12;
  JB->rowvals[26] = 16;
  JB->rowvals[27] = 17;
  JB->rowvals[28] = 6;
  JB->rowvals[29] = 7;
  JB->rowvals[30] = 11;
  JB->rowvals[31] = 14;
  JB->rowvals[32] = 7;
  JB->rowvals[33] = 11;
  JB->rowvals[34] = 7;
  JB->rowvals[35] = 8;
  JB->rowvals[36] = 9;
  JB->rowvals[37] = 17;
  JB->rowvals[38] = 6;
  JB->rowvals[39] = 9;
  JB->rowvals[40] = 16;
  JB->rowvals[41] = 10;
  JB->rowvals[42] = 11;
  JB->rowvals[43] = 17;
  JB->rowvals[44] = 10;
  JB->rowvals[45] = 11;
  JB->rowvals[46] = 17;
  JB->rowvals[47] = 4;
  JB->rowvals[48] = 5;
  JB->rowvals[49] = 10;
  JB->rowvals[50] = 11;
  JB->rowvals[51] = 12;
  JB->rowvals[52] = 15;
  JB->rowvals[53] = 16;
  JB->rowvals[54] = 17;
  JB->rowvals[55] = 9;
  JB->rowvals[56] = 10;
  JB->rowvals[57] = 12;
  JB->rowvals[58] = 13;
  JB->rowvals[59] = 14;
  JB->rowvals[60] = 17;
  JB->rowvals[61] = 5;
  JB->rowvals[62] = 9;
  JB->rowvals[63] = 10;
  JB->rowvals[64] = 13;
  JB->rowvals[65] = 14;
  JB->rowvals[66] = 17;
  JB->rowvals[67] = 8;
  JB->rowvals[68] = 10;
  JB->rowvals[69] = 13;
  JB->rowvals[70] = 15;
  JB->rowvals[71] = 16;
  JB->rowvals[72] = 17;
  JB->rowvals[73] = 8;
  JB->rowvals[74] = 10;
  JB->rowvals[75] = 14;
  JB->rowvals[76] = 15;
  JB->rowvals[77] = 16;
  JB->rowvals[78] = 17;
  JB->rowvals[79] = 7;
  JB->rowvals[80] = 17;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 2;
  JB->colptrs[2] = 4;
  JB->colptrs[3] = 10;
  JB->colptrs[4] = 16;
  JB->colptrs[5] = 22;
  JB->colptrs[6] = 28;
  JB->colptrs[7] = 32;
  JB->colptrs[8] = 34;
  JB->colptrs[9] = 38;
  JB->colptrs[10] = 41;
  JB->colptrs[11] = 44;
  JB->colptrs[12] = 47;
  JB->colptrs[13] = 55;
  JB->colptrs[14] = 61;
  JB->colptrs[15] = 67;
  JB->colptrs[16] = 73;
  JB->colptrs[17] = 79;
  JB->colptrs[18] = 81;
  return(0);
}


 int sx_geneExp_EMRE(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*18);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + x_tmp[3] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] - (1.0*(2.0*k[0]*x_tmp[4] - 1.0*x_tmp[10]))/k[0] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (2.0*k[0]*x_tmp[12] + x_tmp[10])/k[0] - 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + x_tmp[10] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*(k[0]*x_tmp[10] - 1.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[13] + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*x_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[15] + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 1: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] - 1.0*x_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[2] - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = (2.0*k[0]*x_tmp[12] + x_tmp[11])/k[0] + 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - (1.0*(2.0*k[0]*x_tmp[5] - 1.0*x_tmp[11]))/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[11] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] - 1.0*x_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*(k[0]*x_tmp[11] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] + x_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[14] + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] + x_tmp[16] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[16] + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 2: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0] + x_tmp[2];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = (2.0*k[0]*x_tmp[14] + x_tmp[11])/k[0] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[14] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7] + x_tmp[11];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]) + x_tmp[16];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] + x_tmp[12] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[5] + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 3: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0] - 1.0*x_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] - (1.0*(2.0*k[0]*x_tmp[6] - 1.0*x_tmp[7]))/k[0] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7] - 1.0*x_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]) - 1.0*x_tmp[9];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[14] + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 4: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1] + x_tmp[0];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = (2.0*k[0]*x_tmp[9] + x_tmp[7])/k[0] + 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]) + x_tmp[6];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] + x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + x_tmp[14] + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17] + x_tmp[7];

  } break;

  case 5: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1] - 1.0*x_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] - (1.0*(2.0*k[0]*x_tmp[8] - 1.0*x_tmp[17]))/k[0] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]) - 1.0*x_tmp[9];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*x_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*x_tmp[16] + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17] - 1.0*x_tmp[17];

  } break;

  case 6: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = k[0]*x_tmp[15] - 1.0*p[1]*sx_tmp[2] + sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[15] + k[0]*x_tmp[1]*x_tmp[10] + k[0]*x_tmp[3]*x_tmp[17] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*k[0]*x_tmp[15] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*x_tmp[1]*x_tmp[10] - 1.0*k[0]*x_tmp[3]*x_tmp[17] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] - (1.0*(2.0*(pow(k[0],2))*x_tmp[4]*x_tmp[17] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[15] - 1.0*k[0]*x_tmp[10]*x_tmp[17]))/k[0] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[12]*x_tmp[17] + k[0]*x_tmp[10]*x_tmp[17])/k[0] - 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*x_tmp[10]*x_tmp[17] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*x_tmp[10]*x_tmp[17] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*((pow(k[0],2))*x_tmp[10]*x_tmp[17] - 1.0*(pow(k[0],3))*x_tmp[10]*x_tmp[15] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[17] + (pow(k[0],3))*x_tmp[10]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[17]))/(pow(k[0],2)) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*x_tmp[9]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*x_tmp[17] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*x_tmp[9]*x_tmp[10] + k[0]*x_tmp[13]*x_tmp[17] + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[15]*x_tmp[17] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*x_tmp[8]*x_tmp[10] + k[0]*x_tmp[15]*x_tmp[17] + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 7: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 8: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  case 9: {
sxdot_tmp[0] = p[2]*sx_tmp[2] - 1.0*p[3]*sx_tmp[0];
sxdot_tmp[1] = p[4]*sx_tmp[0] - 1.0*p[5]*sx_tmp[1];
sxdot_tmp[2] = sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[2] + k[0]*p[6]*sx_tmp[15] + k[0]*p[6]*sx_tmp[1]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[1] + k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[3] = p[1]*sx_tmp[2] - 1.0*sx_tmp[3]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[15] - 1.0*k[0]*p[6]*sx_tmp[1]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[1] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[3];
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[12] + (p[1]*sx_tmp[11])/k[0] + (sx_tmp[10]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[15] + k[0]*p[6]*x_tmp[17]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[4])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[4] - 1.0*k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] - 2.0*k[0]*p[6]*sx_tmp[15]*x_tmp[10];
sxdot_tmp[5] = (p[1]*sx_tmp[11])/k[0] - 2.0*p[1]*sx_tmp[5] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[17])*sx_tmp[12])/k[0] + (sx_tmp[10]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[16] + k[0]*p[6]*x_tmp[17]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[10])*sx_tmp[17])/k[0] + 2.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[14] - 2.0*p[3]*sx_tmp[6] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[11])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[7];
sxdot_tmp[8] = 2.0*p[4]*sx_tmp[9] - 2.0*p[5]*sx_tmp[8] + (p[4]*sx_tmp[7])/k[0] + (p[5]*sx_tmp[17])/k[0];
sxdot_tmp[9] = p[4]*sx_tmp[6] + p[2]*sx_tmp[16] - 1.0*sx_tmp[9]*(p[3] + p[5]);
sxdot_tmp[10] = p[1]*sx_tmp[11] - 1.0*sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[11] = sx_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[17]) - 1.0*p[1]*sx_tmp[11] + k[0]*p[6]*sx_tmp[17]*x_tmp[10];
sxdot_tmp[12] = p[1]*sx_tmp[5] + sx_tmp[4]*(p[0] + k[0]*p[6]*x_tmp[17]) - (1.0*sx_tmp[10]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[15] + (pow(k[0],2))*p[6]*x_tmp[17] + (pow(k[0],3))*p[6]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[6]*x_tmp[10] - 1.0*(pow(k[0],3))*p[6]*x_tmp[4] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[17]))*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[15]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[16]*x_tmp[10];
sxdot_tmp[13] = p[2]*sx_tmp[12] - 1.0*sx_tmp[13]*(p[0] + p[3] + k[0]*p[6]*x_tmp[17]) + p[1]*sx_tmp[14] - 1.0*k[0]*p[6]*sx_tmp[9]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] - 1.0*sx_tmp[14]*(p[1] + p[3]) + sx_tmp[13]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[9]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[9] + k[0]*p[6]*sx_tmp[17]*x_tmp[13];
sxdot_tmp[15] = p[1]*sx_tmp[16] - 1.0*sx_tmp[15]*(p[0] + p[5] + k[0]*p[6]*x_tmp[17]) + p[4]*sx_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[8]*x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[10]*x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[16] = p[4]*sx_tmp[14] - 1.0*sx_tmp[16]*(p[1] + p[5]) + sx_tmp[15]*(p[0] + k[0]*p[6]*x_tmp[17]) + k[0]*p[6]*sx_tmp[8]*x_tmp[10] + k[0]*p[6]*sx_tmp[10]*x_tmp[8] + k[0]*p[6]*sx_tmp[17]*x_tmp[15];
sxdot_tmp[17] = p[4]*sx_tmp[7] - 1.0*p[5]*sx_tmp[17];

  } break;

  }
 for (ix=0; ix<18; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_geneExp_EMRE(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*18);
  switch (ip) {
  case 9: {
sx0_tmp[7] = -(1.0*(k[3] - 1.0))/k[0];

  } break;

  }

  return;
}


void y_geneExp_EMRE(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*3] + x[it+nt*10];
y[it+nt*1] = x[it+nt*2] + x[it+nt*11];
y[it+nt*2] = x[it+nt*0] + x[it+nt*7];
y[it+nt*3] = x[it+nt*1] + x[it+nt*17];
y[it+nt*4] = x[it+nt*4];
y[it+nt*5] = x[it+nt*12];
y[it+nt*6] = x[it+nt*5];
y[it+nt*7] = x[it+nt*13];
y[it+nt*8] = x[it+nt*14];
y[it+nt*9] = x[it+nt*6];
y[it+nt*10] = x[it+nt*15];
y[it+nt*11] = x[it+nt*16];
y[it+nt*12] = x[it+nt*9];
y[it+nt*13] = x[it+nt*8];
y[it+nt*14] = (x[it+nt*1] + x[it+nt*17])*p[7] + p[8];
y[it+nt*15] = (pow(p[7],2))*x[it+nt*8];
    
    return;
}


void dydp_geneExp_EMRE(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  memset(dydp,0,sizeof(double)*16*np*nt);
  switch (plist[ip]) {
  case 7: {
dydp[it+nt*(14+ip*ny)] = x[it+nt*1] + x[it+nt*17];
dydp[it+nt*(15+ip*ny)] = 2.0*p[7]*x[it+nt*8];

  } break;

  case 8: {
dydp[it+nt*(14+ip*ny)] = 1.0;

  } break;

  }
  }
  
  return;
}


void dydx_geneExp_EMRE(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*288);
dydx[2] = 1.0;
dydx[19] = 1.0;
dydx[30] = p[7];
dydx[33] = 1.0;
dydx[48] = 1.0;
dydx[68] = 1.0;
dydx[86] = 1.0;
dydx[105] = 1.0;
dydx[114] = 1.0;
dydx[141] = 1.0;
dydx[143] = (pow(p[7],2));
dydx[156] = 1.0;
dydx[160] = 1.0;
dydx[177] = 1.0;
dydx[197] = 1.0;
dydx[215] = 1.0;
dydx[232] = 1.0;
dydx[250] = 1.0;
dydx[267] = 1.0;
dydx[275] = 1.0;
dydx[286] = p[7];
  
  return;
}


void sy_geneExp_EMRE(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)] + x[it+nt*1] + x[it+nt*17];
sy[it+nt*(15+np*ny)] = 2.0*p[7]*x[it+nt*8] + (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 8: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)] + 1.0;
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  case 9: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(3+np*nx)] + sx[it+nt*(10+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(2+np*nx)] + sx[it+nt*(11+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(0+np*nx)] + sx[it+nt*(7+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(1+np*nx)] + sx[it+nt*(17+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(14+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(15+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(16+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(1+np*nx)] + p[7]*sx[it+nt*(17+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(8+np*nx)];

  } break;

  }
  
  return;
}
int root_geneExp_EMRE(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_geneExp_EMRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_geneExp_EMRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_geneExp_EMRE(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_geneExp_EMRE(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
