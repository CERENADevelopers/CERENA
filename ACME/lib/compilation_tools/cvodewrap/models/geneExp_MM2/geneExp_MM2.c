#include "geneExp_MM2.h"
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


 int xdot_geneExp_MM2(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*14);
xdot_tmp[0] = (2.0*(pow(k[0],2))*p[1]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[6]*x_tmp[13] + k[0]*p[0]*x_tmp[5] + k[0]*p[1]*x_tmp[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[0]*x_tmp[7] + (pow(k[0],2))*p[6]*x_tmp[5]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[1] = (2.0*(pow(k[0],2))*p[2]*x_tmp[11] - 2.0*(pow(k[0],2))*p[3]*x_tmp[1] + k[0]*p[3]*x_tmp[2] + k[0]*p[2]*x_tmp[6])/(pow(k[0],2));
xdot_tmp[2] = -(1.0*(k[0]*p[3]*x_tmp[2] - 1.0*k[0]*p[2]*x_tmp[6]))/k[0];
xdot_tmp[3] = (2.0*(pow(k[0],2))*p[4]*x_tmp[4] - 2.0*(pow(k[0],2))*p[5]*x_tmp[3] + k[0]*p[4]*x_tmp[2] + k[0]*p[5]*x_tmp[7])/(pow(k[0],2));
xdot_tmp[4] = ((pow(k[0],2))*p[4]*x_tmp[1] - 1.0*(pow(k[0],2))*p[3]*x_tmp[4] - 1.0*(pow(k[0],2))*p[5]*x_tmp[4] + (pow(k[0],2))*p[2]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[5] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[13] + k[0]*p[0]*x_tmp[5] - 1.0*k[0]*p[1]*x_tmp[6] + (pow(k[0],2))*p[6]*x_tmp[5]*x_tmp[7]))/k[0];
xdot_tmp[6] = ((pow(k[0],2))*p[6]*x_tmp[13] + k[0]*p[0]*x_tmp[5] - 1.0*k[0]*p[1]*x_tmp[6] + (pow(k[0],2))*p[6]*x_tmp[5]*x_tmp[7])/k[0];
xdot_tmp[7] = (k[0]*p[4]*x_tmp[2] - 1.0*k[0]*p[5]*x_tmp[7])/k[0];
xdot_tmp[8] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[1]*x_tmp[8] - 1.0*(pow(k[0],2))*p[1]*x_tmp[10] + (pow(k[0],2))*p[6]*x_tmp[13] + k[0]*p[0]*x_tmp[5] + k[0]*p[1]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0]*x_tmp[7] + (pow(k[0],2))*p[6]*x_tmp[5]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[8] + (pow(k[0],3))*p[6]*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5]*x_tmp[13]))/(pow(k[0],2));
xdot_tmp[9] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[2]*x_tmp[8] - 1.0*(pow(k[0],2))*p[1]*x_tmp[11] + (pow(k[0],2))*p[3]*x_tmp[9] + (pow(k[0],3))*p[6]*x_tmp[4]*x_tmp[5] + (pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[9]))/(pow(k[0],2));
xdot_tmp[10] = (2.0*(pow(k[0],2))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[1]*x_tmp[10] + (pow(k[0],2))*p[6]*x_tmp[13] + k[0]*p[0]*x_tmp[5] + k[0]*p[1]*x_tmp[6] + (pow(k[0],2))*p[6]*x_tmp[5]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[8] + 2.0*(pow(k[0],3))*p[6]*x_tmp[5]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[11] = ((pow(k[0],2))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[1]*x_tmp[11] + (pow(k[0],2))*p[2]*x_tmp[10] - 1.0*(pow(k[0],2))*p[3]*x_tmp[11] + (pow(k[0],3))*p[6]*x_tmp[4]*x_tmp[5] + (pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[9])/(pow(k[0],2));
xdot_tmp[12] = ((pow(k[0],2))*p[0]*x_tmp[13] - 1.0*(pow(k[0],2))*p[1]*x_tmp[12] + (pow(k[0],2))*p[4]*x_tmp[11] - 1.0*(pow(k[0],2))*p[5]*x_tmp[12] + (pow(k[0],3))*p[6]*x_tmp[3]*x_tmp[5] + (pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[13] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[13] - 1.0*(pow(k[0],2))*p[1]*x_tmp[12] - 1.0*(pow(k[0],2))*p[4]*x_tmp[9] + (pow(k[0],2))*p[5]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[3]*x_tmp[5] + (pow(k[0],3))*p[6]*x_tmp[7]*x_tmp[13]))/(pow(k[0],2));

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_geneExp_MM2(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*28);

  for (ixB=0; ixB<28; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_geneExp_MM2(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
  memset(qBdot_tmp,0,sizeof(double)*2*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_geneExp_MM2(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*14);
x0_tmp[0] = (k[5]*k[19])/(pow(k[0],2));
x0_tmp[1] = (k[12]*k[26])/(pow(k[0],2));
x0_tmp[2] = -(1.0*(p[9]*(k[3] - 1.0) - 1.0*k[3]*k[17]))/k[0];
x0_tmp[3] = (k[14]*k[28])/(pow(k[0],2));
x0_tmp[4] = (k[13]*k[27])/(pow(k[0],2));
x0_tmp[5] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[6] = (k[2]*k[16])/k[0];
x0_tmp[7] = (k[4]*k[18])/k[0];
x0_tmp[8] = (k[6]*k[20])/(pow(k[0],2));
x0_tmp[9] = (k[7]*k[21])/(pow(k[0],2));
x0_tmp[10] = (k[9]*k[23])/(pow(k[0],2));
x0_tmp[11] = (k[10]*k[24])/(pow(k[0],2));
x0_tmp[12] = (k[11]*k[25])/(pow(k[0],2));
x0_tmp[13] = (k[8]*k[22])/(pow(k[0],2));
  
  
  return;
}


 int Jv_geneExp_MM2(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*14);
Jv_tmp[0] = 2.0*p[1]*v_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*v_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*v_tmp[7])/(pow(k[0],2)) + (p[1]*v_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*v_tmp[0])/(pow(k[0],2)) + (v_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
Jv_tmp[1] = 2.0*p[2]*v_tmp[11] - 2.0*p[3]*v_tmp[1] + (p[3]*v_tmp[2])/k[0] + (p[2]*v_tmp[6])/k[0];
Jv_tmp[2] = p[2]*v_tmp[6] - 1.0*p[3]*v_tmp[2];
Jv_tmp[3] = 2.0*p[4]*v_tmp[4] - 2.0*p[5]*v_tmp[3] + (p[4]*v_tmp[2])/k[0] + (p[5]*v_tmp[7])/k[0];
Jv_tmp[4] = p[4]*v_tmp[1] + p[2]*v_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*v_tmp[4])/(pow(k[0],2));
Jv_tmp[5] = p[1]*v_tmp[6] - 1.0*k[0]*p[6]*v_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*v_tmp[5])/k[0] - 1.0*k[0]*p[6]*v_tmp[7]*x_tmp[5];
Jv_tmp[6] = k[0]*p[6]*v_tmp[13] - 1.0*p[1]*v_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*v_tmp[5])/k[0] + k[0]*p[6]*v_tmp[7]*x_tmp[5];
Jv_tmp[7] = p[4]*v_tmp[2] - 1.0*p[5]*v_tmp[7];
Jv_tmp[8] = p[1]*v_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*v_tmp[13])/(pow(k[0],2)) - (1.0*v_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*v_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*v_tmp[6])/k[0] + (v_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*v_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*v_tmp[12]*x_tmp[5];
Jv_tmp[9] = p[2]*v_tmp[8] + p[1]*v_tmp[11] - (1.0*v_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*v_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*v_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*v_tmp[7]*x_tmp[9];
Jv_tmp[10] = p[6]*v_tmp[13] - 2.0*p[1]*v_tmp[10] + (p[1]*v_tmp[6])/k[0] + (v_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (v_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*v_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*v_tmp[12]*x_tmp[5];
Jv_tmp[11] = p[2]*v_tmp[10] + (v_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*v_tmp[4]*x_tmp[5] + k[0]*p[6]*v_tmp[5]*x_tmp[4] + k[0]*p[6]*v_tmp[7]*x_tmp[9];
Jv_tmp[12] = p[4]*v_tmp[11] + (v_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*v_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*v_tmp[3]*x_tmp[5] + k[0]*p[6]*v_tmp[5]*x_tmp[3] + k[0]*p[6]*v_tmp[7]*x_tmp[13];
Jv_tmp[13] = p[1]*v_tmp[12] + p[4]*v_tmp[9] - (1.0*v_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*v_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*v_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*v_tmp[7]*x_tmp[13];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_geneExp_MM2(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*14);
JvB_tmp[0] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*vB_tmp[0])/(pow(k[0],2)) - (1.0*vB_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
JvB_tmp[1] = 2.0*p[3]*vB_tmp[1] - 1.0*p[4]*vB_tmp[4];
JvB_tmp[2] = p[3]*vB_tmp[2] - 1.0*p[4]*vB_tmp[7] - (1.0*p[3]*vB_tmp[1])/k[0] - (1.0*p[4]*vB_tmp[3])/k[0];
JvB_tmp[3] = 2.0*p[5]*vB_tmp[3] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[5] + k[0]*p[6]*vB_tmp[13]*x_tmp[5];
JvB_tmp[4] = (((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*vB_tmp[4])/(pow(k[0],2)) - 2.0*p[4]*vB_tmp[3] + k[0]*p[6]*vB_tmp[9]*x_tmp[5] - 1.0*k[0]*p[6]*vB_tmp[11]*x_tmp[5];
JvB_tmp[5] = (vB_tmp[8]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[10]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*vB_tmp[0]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*vB_tmp[5])/k[0] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*vB_tmp[6])/k[0] + k[0]*p[6]*vB_tmp[9]*x_tmp[4] - 1.0*k[0]*p[6]*vB_tmp[11]*x_tmp[4] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[3] + k[0]*p[6]*vB_tmp[13]*x_tmp[3];
JvB_tmp[6] = p[1]*vB_tmp[6] - 1.0*p[1]*vB_tmp[5] - 1.0*p[2]*vB_tmp[2] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[2]*vB_tmp[1])/k[0] + (p[1]*vB_tmp[8])/k[0] - (1.0*p[1]*vB_tmp[10])/k[0];
JvB_tmp[7] = p[5]*vB_tmp[7] + (vB_tmp[8]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*vB_tmp[0])/(pow(k[0],2)) - (1.0*p[5]*vB_tmp[3])/k[0] - (1.0*vB_tmp[10]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[6]*vB_tmp[5]*x_tmp[5] - 1.0*k[0]*p[6]*vB_tmp[6]*x_tmp[5] + k[0]*p[6]*vB_tmp[9]*x_tmp[9] - 1.0*k[0]*p[6]*vB_tmp[11]*x_tmp[9] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[13] + k[0]*p[6]*vB_tmp[13]*x_tmp[13];
JvB_tmp[8] = (vB_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[9] - 2.0*p[1]*vB_tmp[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*vB_tmp[10])/(pow(k[0],2));
JvB_tmp[9] = (vB_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*vB_tmp[11]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*p[4]*vB_tmp[13];
JvB_tmp[10] = 2.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[8] - 1.0*p[2]*vB_tmp[11];
JvB_tmp[11] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[11])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[9] - 1.0*p[4]*vB_tmp[12] - 2.0*p[2]*vB_tmp[1];
JvB_tmp[12] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*vB_tmp[12])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[13] - 1.0*p[2]*vB_tmp[4] + k[0]*p[6]*vB_tmp[8]*x_tmp[5] - 2.0*k[0]*p[6]*vB_tmp[10]*x_tmp[5];
JvB_tmp[13] = (((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*vB_tmp[8])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*vB_tmp[0])/(pow(k[0],2)) - 1.0*p[6]*vB_tmp[10] - (1.0*vB_tmp[12]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) + (vB_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) + k[0]*p[6]*vB_tmp[5] - 1.0*k[0]*p[6]*vB_tmp[6];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_geneExp_MM2(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_geneExp_MM2(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_geneExp_MM2(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
J->data[0] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[8] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[15] = -2.0*p[3];
J->data[18] = p[4];
J->data[29] = p[3]/k[0];
J->data[30] = -1.0*p[3];
J->data[31] = p[4]/k[0];
J->data[35] = p[4];
J->data[45] = -2.0*p[5];
J->data[54] = k[0]*p[6]*x_tmp[5];
J->data[55] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[59] = 2.0*p[4];
J->data[60] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5]))/(pow(k[0],2));
J->data[65] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[67] = k[0]*p[6]*x_tmp[5];
J->data[70] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13])/(pow(k[0],2));
J->data[75] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7]))/k[0];
J->data[76] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])/k[0];
J->data[78] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
J->data[79] = -1.0*k[0]*p[6]*x_tmp[4];
J->data[80] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12])/(pow(k[0],2));
J->data[81] = k[0]*p[6]*x_tmp[4];
J->data[82] = k[0]*p[6]*x_tmp[3];
J->data[83] = -1.0*k[0]*p[6]*x_tmp[3];
J->data[84] = p[1]/k[0];
J->data[85] = p[2]/k[0];
J->data[86] = p[2];
J->data[89] = p[1];
J->data[90] = -1.0*p[1];
J->data[92] = -(1.0*p[1])/k[0];
J->data[94] = p[1]/k[0];
J->data[98] = -(1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5]))/(pow(k[0],2));
J->data[101] = p[5]/k[0];
J->data[103] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[104] = k[0]*p[6]*x_tmp[5];
J->data[105] = -1.0*p[5];
J->data[106] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2));
J->data[107] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[108] = ((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8])/(pow(k[0],2));
J->data[109] = k[0]*p[6]*x_tmp[9];
J->data[110] = k[0]*p[6]*x_tmp[13];
J->data[111] = -1.0*k[0]*p[6]*x_tmp[13];
J->data[112] = 2.0*p[1];
J->data[120] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[121] = p[2];
J->data[122] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[135] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[137] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[139] = p[4];
J->data[148] = p[1];
J->data[150] = -2.0*p[1];
J->data[151] = p[2];
J->data[155] = 2.0*p[2];
J->data[163] = p[1];
J->data[165] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[166] = p[4];
J->data[172] = p[2];
J->data[176] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[178] = 2.0*k[0]*p[6]*x_tmp[5];
J->data[180] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5]))/(pow(k[0],2));
J->data[181] = p[1];
J->data[182] = ((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])/(pow(k[0],2));
J->data[187] = -1.0*k[0]*p[6];
J->data[188] = k[0]*p[6];
J->data[190] = -(1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5]))/(pow(k[0],2));
J->data[192] = p[6];
J->data[194] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[195] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_geneExp_MM2(realtype t, N_Vector x,
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
  J->rowvals[1] = 8;
  J->rowvals[2] = 1;
  J->rowvals[3] = 4;
  J->rowvals[4] = 1;
  J->rowvals[5] = 2;
  J->rowvals[6] = 3;
  J->rowvals[7] = 7;
  J->rowvals[8] = 3;
  J->rowvals[9] = 12;
  J->rowvals[10] = 13;
  J->rowvals[11] = 3;
  J->rowvals[12] = 4;
  J->rowvals[13] = 9;
  J->rowvals[14] = 11;
  J->rowvals[15] = 0;
  J->rowvals[16] = 5;
  J->rowvals[17] = 6;
  J->rowvals[18] = 8;
  J->rowvals[19] = 9;
  J->rowvals[20] = 10;
  J->rowvals[21] = 11;
  J->rowvals[22] = 12;
  J->rowvals[23] = 13;
  J->rowvals[24] = 0;
  J->rowvals[25] = 1;
  J->rowvals[26] = 2;
  J->rowvals[27] = 5;
  J->rowvals[28] = 6;
  J->rowvals[29] = 8;
  J->rowvals[30] = 10;
  J->rowvals[31] = 0;
  J->rowvals[32] = 3;
  J->rowvals[33] = 5;
  J->rowvals[34] = 6;
  J->rowvals[35] = 7;
  J->rowvals[36] = 8;
  J->rowvals[37] = 9;
  J->rowvals[38] = 10;
  J->rowvals[39] = 11;
  J->rowvals[40] = 12;
  J->rowvals[41] = 13;
  J->rowvals[42] = 0;
  J->rowvals[43] = 8;
  J->rowvals[44] = 9;
  J->rowvals[45] = 10;
  J->rowvals[46] = 9;
  J->rowvals[47] = 11;
  J->rowvals[48] = 13;
  J->rowvals[49] = 8;
  J->rowvals[50] = 10;
  J->rowvals[51] = 11;
  J->rowvals[52] = 1;
  J->rowvals[53] = 9;
  J->rowvals[54] = 11;
  J->rowvals[55] = 12;
  J->rowvals[56] = 4;
  J->rowvals[57] = 8;
  J->rowvals[58] = 10;
  J->rowvals[59] = 12;
  J->rowvals[60] = 13;
  J->rowvals[61] = 0;
  J->rowvals[62] = 5;
  J->rowvals[63] = 6;
  J->rowvals[64] = 8;
  J->rowvals[65] = 10;
  J->rowvals[66] = 12;
  J->rowvals[67] = 13;
  J->colptrs[0] = 0;
  J->colptrs[1] = 2;
  J->colptrs[2] = 4;
  J->colptrs[3] = 8;
  J->colptrs[4] = 11;
  J->colptrs[5] = 15;
  J->colptrs[6] = 24;
  J->colptrs[7] = 31;
  J->colptrs[8] = 42;
  J->colptrs[9] = 46;
  J->colptrs[10] = 49;
  J->colptrs[11] = 52;
  J->colptrs[12] = 56;
  J->colptrs[13] = 61;
  J->colptrs[14] = 68;
J->data[0] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[1] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[2] = -2.0*p[3];
J->data[3] = p[4];
J->data[4] = p[3]/k[0];
J->data[5] = -1.0*p[3];
J->data[6] = p[4]/k[0];
J->data[7] = p[4];
J->data[8] = -2.0*p[5];
J->data[9] = k[0]*p[6]*x_tmp[5];
J->data[10] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[11] = 2.0*p[4];
J->data[12] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5]))/(pow(k[0],2));
J->data[13] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[14] = k[0]*p[6]*x_tmp[5];
J->data[15] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13])/(pow(k[0],2));
J->data[16] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7]))/k[0];
J->data[17] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])/k[0];
J->data[18] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
J->data[19] = -1.0*k[0]*p[6]*x_tmp[4];
J->data[20] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12])/(pow(k[0],2));
J->data[21] = k[0]*p[6]*x_tmp[4];
J->data[22] = k[0]*p[6]*x_tmp[3];
J->data[23] = -1.0*k[0]*p[6]*x_tmp[3];
J->data[24] = p[1]/k[0];
J->data[25] = p[2]/k[0];
J->data[26] = p[2];
J->data[27] = p[1];
J->data[28] = -1.0*p[1];
J->data[29] = -(1.0*p[1])/k[0];
J->data[30] = p[1]/k[0];
J->data[31] = -(1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5]))/(pow(k[0],2));
J->data[32] = p[5]/k[0];
J->data[33] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[34] = k[0]*p[6]*x_tmp[5];
J->data[35] = -1.0*p[5];
J->data[36] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2));
J->data[37] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[38] = ((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8])/(pow(k[0],2));
J->data[39] = k[0]*p[6]*x_tmp[9];
J->data[40] = k[0]*p[6]*x_tmp[13];
J->data[41] = -1.0*k[0]*p[6]*x_tmp[13];
J->data[42] = 2.0*p[1];
J->data[43] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[44] = p[2];
J->data[45] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[46] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
J->data[47] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[48] = p[4];
J->data[49] = p[1];
J->data[50] = -2.0*p[1];
J->data[51] = p[2];
J->data[52] = 2.0*p[2];
J->data[53] = p[1];
J->data[54] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[55] = p[4];
J->data[56] = p[2];
J->data[57] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[58] = 2.0*k[0]*p[6]*x_tmp[5];
J->data[59] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5]))/(pow(k[0],2));
J->data[60] = p[1];
J->data[61] = ((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])/(pow(k[0],2));
J->data[62] = -1.0*k[0]*p[6];
J->data[63] = k[0]*p[6];
J->data[64] = -(1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5]))/(pow(k[0],2));
J->data[65] = p[6];
J->data[66] = ((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
J->data[67] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
  return(0);
}


 int JBBand_geneExp_MM2(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_geneExp_MM2(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_geneExp_MM2(long int N, realtype t, N_Vector x,
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
JB->data[0] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
JB->data[5] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
JB->data[6] = -(1.0*p[1])/k[0];
JB->data[7] = (2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])/(pow(k[0],2));
JB->data[8] = -2.0*p[1];
JB->data[13] = -(1.0*((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5]))/(pow(k[0],2));
JB->data[15] = 2.0*p[3];
JB->data[16] = -(1.0*p[3])/k[0];
JB->data[20] = -(1.0*p[2])/k[0];
JB->data[25] = -2.0*p[2];
JB->data[30] = p[3];
JB->data[34] = -1.0*p[2];
JB->data[44] = -(1.0*p[4])/k[0];
JB->data[45] = 2.0*p[5];
JB->data[46] = -2.0*p[4];
JB->data[49] = -(1.0*p[5])/k[0];
JB->data[57] = -1.0*p[4];
JB->data[60] = ((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])/(pow(k[0],2));
JB->data[68] = -1.0*p[2];
JB->data[75] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])/k[0];
JB->data[76] = -1.0*p[1];
JB->data[77] = k[0]*p[6]*x_tmp[5];
JB->data[83] = k[0]*p[6];
JB->data[89] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7]))/k[0];
JB->data[90] = p[1];
JB->data[91] = -1.0*k[0]*p[6]*x_tmp[5];
JB->data[97] = -1.0*k[0]*p[6];
JB->data[100] = -1.0*p[4];
JB->data[105] = p[5];
JB->data[112] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
JB->data[117] = (k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13])/(pow(k[0],2));
JB->data[118] = p[1]/k[0];
JB->data[119] = ((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8])/(pow(k[0],2));
JB->data[120] = ((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
JB->data[122] = -1.0*p[1];
JB->data[124] = k[0]*p[6]*x_tmp[5];
JB->data[125] = ((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])/(pow(k[0],2));
JB->data[130] = k[0]*p[6]*x_tmp[5];
JB->data[131] = k[0]*p[6]*x_tmp[4];
JB->data[133] = k[0]*p[6]*x_tmp[9];
JB->data[134] = -1.0*p[2];
JB->data[135] = ((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));
JB->data[137] = -1.0*p[1];
JB->data[145] = -(1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2));
JB->data[146] = -(1.0*p[1])/k[0];
JB->data[147] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2));
JB->data[148] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
JB->data[150] = 2.0*p[1];
JB->data[152] = -2.0*k[0]*p[6]*x_tmp[5];
JB->data[153] = -1.0*p[6];
JB->data[158] = -1.0*k[0]*p[6]*x_tmp[5];
JB->data[159] = -1.0*k[0]*p[6]*x_tmp[4];
JB->data[161] = -1.0*k[0]*p[6]*x_tmp[9];
JB->data[163] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
JB->data[164] = -1.0*p[2];
JB->data[165] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[171] = -1.0*k[0]*p[6]*x_tmp[5];
JB->data[173] = -1.0*k[0]*p[6]*x_tmp[3];
JB->data[175] = -1.0*k[0]*p[6]*x_tmp[13];
JB->data[179] = -1.0*p[4];
JB->data[180] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])/(pow(k[0],2));
JB->data[181] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2));
JB->data[185] = k[0]*p[6]*x_tmp[5];
JB->data[187] = k[0]*p[6]*x_tmp[3];
JB->data[189] = k[0]*p[6]*x_tmp[13];
JB->data[191] = -1.0*p[4];
JB->data[194] = -1.0*p[1];
JB->data[195] = ((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7])/(pow(k[0],2));

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_geneExp_MM2(realtype t, N_Vector x,
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
  JB->rowvals[1] = 5;
  JB->rowvals[2] = 6;
  JB->rowvals[3] = 7;
  JB->rowvals[4] = 8;
  JB->rowvals[5] = 13;
  JB->rowvals[6] = 1;
  JB->rowvals[7] = 2;
  JB->rowvals[8] = 6;
  JB->rowvals[9] = 11;
  JB->rowvals[10] = 2;
  JB->rowvals[11] = 6;
  JB->rowvals[12] = 2;
  JB->rowvals[13] = 3;
  JB->rowvals[14] = 4;
  JB->rowvals[15] = 7;
  JB->rowvals[16] = 1;
  JB->rowvals[17] = 4;
  JB->rowvals[18] = 12;
  JB->rowvals[19] = 5;
  JB->rowvals[20] = 6;
  JB->rowvals[21] = 7;
  JB->rowvals[22] = 13;
  JB->rowvals[23] = 5;
  JB->rowvals[24] = 6;
  JB->rowvals[25] = 7;
  JB->rowvals[26] = 13;
  JB->rowvals[27] = 2;
  JB->rowvals[28] = 7;
  JB->rowvals[29] = 0;
  JB->rowvals[30] = 5;
  JB->rowvals[31] = 6;
  JB->rowvals[32] = 7;
  JB->rowvals[33] = 8;
  JB->rowvals[34] = 10;
  JB->rowvals[35] = 12;
  JB->rowvals[36] = 13;
  JB->rowvals[37] = 4;
  JB->rowvals[38] = 5;
  JB->rowvals[39] = 7;
  JB->rowvals[40] = 8;
  JB->rowvals[41] = 9;
  JB->rowvals[42] = 11;
  JB->rowvals[43] = 5;
  JB->rowvals[44] = 6;
  JB->rowvals[45] = 7;
  JB->rowvals[46] = 8;
  JB->rowvals[47] = 10;
  JB->rowvals[48] = 12;
  JB->rowvals[49] = 13;
  JB->rowvals[50] = 4;
  JB->rowvals[51] = 5;
  JB->rowvals[52] = 7;
  JB->rowvals[53] = 9;
  JB->rowvals[54] = 10;
  JB->rowvals[55] = 11;
  JB->rowvals[56] = 3;
  JB->rowvals[57] = 5;
  JB->rowvals[58] = 7;
  JB->rowvals[59] = 11;
  JB->rowvals[60] = 12;
  JB->rowvals[61] = 13;
  JB->rowvals[62] = 3;
  JB->rowvals[63] = 5;
  JB->rowvals[64] = 7;
  JB->rowvals[65] = 9;
  JB->rowvals[66] = 12;
  JB->rowvals[67] = 13;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 6;
  JB->colptrs[2] = 10;
  JB->colptrs[3] = 12;
  JB->colptrs[4] = 16;
  JB->colptrs[5] = 19;
  JB->colptrs[6] = 23;
  JB->colptrs[7] = 27;
  JB->colptrs[8] = 29;
  JB->colptrs[9] = 37;
  JB->colptrs[10] = 43;
  JB->colptrs[11] = 50;
  JB->colptrs[12] = 56;
  JB->colptrs[13] = 62;
  JB->colptrs[14] = 68;
  return(0);
}


 int sx_geneExp_MM2(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*14);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0])/(pow(k[0],2)) + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = x_tmp[5] - 1.0*p[1]*sx_tmp[6] + k[0]*p[6]*sx_tmp[13] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*(k[0]*x_tmp[5] - 1.0*(pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - 1.0*x_tmp[9] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[8])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + x_tmp[9] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + x_tmp[13] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - 1.0*x_tmp[13] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 1: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (k[0]*x_tmp[6] + 2.0*(pow(k[0],2))*x_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] + x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*x_tmp[6] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*(k[0]*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] + x_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (k[0]*x_tmp[6] - 2.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] - 1.0*x_tmp[11] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] - 1.0*x_tmp[12] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] + x_tmp[12] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 2: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (k[0]*x_tmp[6] + 2.0*(pow(k[0],2))*x_tmp[11])/(pow(k[0],2)) + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2] + x_tmp[6];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] + x_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] + x_tmp[8] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + x_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 3: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (k[0]*x_tmp[2] - 2.0*(pow(k[0],2))*x_tmp[1])/(pow(k[0],2)) + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2] - 1.0*x_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - 1.0*x_tmp[4] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - 1.0*x_tmp[9] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] - 1.0*x_tmp[11] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 4: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (k[0]*x_tmp[2] + 2.0*(pow(k[0],2))*x_tmp[4])/(pow(k[0],2)) + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] + x_tmp[1] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7] + x_tmp[2];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + x_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] + x_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 5: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - 1.0*x_tmp[4] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7] - 1.0*x_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] - 1.0*x_tmp[12] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - 1.0*x_tmp[13] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 6: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + ((pow(k[0],2))*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[5]*x_tmp[7] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[13])/(pow(k[0],2)) + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - (1.0*((pow(k[0],2))*x_tmp[13] + (pow(k[0],2))*x_tmp[5]*x_tmp[7]))/k[0] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = ((pow(k[0],2))*x_tmp[13] + (pow(k[0],2))*x_tmp[5]*x_tmp[7])/k[0] - 1.0*p[1]*sx_tmp[6] + k[0]*p[6]*sx_tmp[13] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[5]*x_tmp[7] + (pow(k[0],3))*x_tmp[7]*x_tmp[8] + (pow(k[0],3))*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*((pow(k[0],3))*x_tmp[4]*x_tmp[5] + (pow(k[0],3))*x_tmp[7]*x_tmp[9]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + ((pow(k[0],2))*x_tmp[13] + (pow(k[0],2))*x_tmp[5]*x_tmp[7] + 2.0*(pow(k[0],3))*x_tmp[7]*x_tmp[8] + 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + ((pow(k[0],3))*x_tmp[4]*x_tmp[5] + (pow(k[0],3))*x_tmp[7]*x_tmp[9])/(pow(k[0],2)) + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + ((pow(k[0],3))*x_tmp[3]*x_tmp[5] + (pow(k[0],3))*x_tmp[7]*x_tmp[13])/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*((pow(k[0],3))*x_tmp[3]*x_tmp[5] + (pow(k[0],3))*x_tmp[7]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 7: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 8: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  case 9: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (((pow(k[0],2))*p[6] - 2.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[6]*x_tmp[0] - 1.0*(pow(k[0],2))*p[6]*x_tmp[5])*sx_tmp[7])/(pow(k[0],2)) + (p[1]*sx_tmp[6])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[0])/(pow(k[0],2)) + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] - 2.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[2]*sx_tmp[11] - 2.0*p[3]*sx_tmp[1] + (p[3]*sx_tmp[2])/k[0] + (p[2]*sx_tmp[6])/k[0];
sxdot_tmp[2] = p[2]*sx_tmp[6] - 1.0*p[3]*sx_tmp[2];
sxdot_tmp[3] = 2.0*p[4]*sx_tmp[4] - 2.0*p[5]*sx_tmp[3] + (p[4]*sx_tmp[2])/k[0] + (p[5]*sx_tmp[7])/k[0];
sxdot_tmp[4] = p[4]*sx_tmp[1] + p[2]*sx_tmp[12] - (1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = p[1]*sx_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13] - (1.0*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[6] = k[0]*p[6]*sx_tmp[13] - 1.0*p[1]*sx_tmp[6] + ((k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7])*sx_tmp[5])/k[0] + k[0]*p[6]*sx_tmp[7]*x_tmp[5];
sxdot_tmp[7] = p[4]*sx_tmp[2] - 1.0*p[5]*sx_tmp[7];
sxdot_tmp[8] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*p[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[5])*sx_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + (pow(k[0],3))*p[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[6]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[6])/k[0] + (sx_tmp[0]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[9] = p[2]*sx_tmp[8] + p[1]*sx_tmp[11] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[10] = p[6]*sx_tmp[13] - 2.0*p[1]*sx_tmp[10] + (p[1]*sx_tmp[6])/k[0] + (sx_tmp[5]*(k[0]*p[0] + (pow(k[0],2))*p[6]*x_tmp[7] + 2.0*(pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[7]*((pow(k[0],2))*p[6]*x_tmp[5] + 2.0*(pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[6]*x_tmp[7])*sx_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[5];
sxdot_tmp[11] = p[2]*sx_tmp[10] + (sx_tmp[9]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[11])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[4]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[4] + k[0]*p[6]*sx_tmp[7]*x_tmp[9];
sxdot_tmp[12] = p[4]*sx_tmp[11] + (sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5])*sx_tmp[12])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[3]*x_tmp[5] + k[0]*p[6]*sx_tmp[5]*x_tmp[3] + k[0]*p[6]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[13] = p[1]*sx_tmp[12] + p[4]*sx_tmp[9] - (1.0*sx_tmp[13]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[6]*x_tmp[7]))/(pow(k[0],2)) - 1.0*k[0]*p[6]*sx_tmp[3]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[3] - 1.0*k[0]*p[6]*sx_tmp[7]*x_tmp[13];

  } break;

  }
 for (ix=0; ix<14; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_geneExp_MM2(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*14);
  switch (ip) {
  case 9: {
sx0_tmp[2] = -(1.0*(k[3] - 1.0))/k[0];

  } break;

  }

  return;
}


void y_geneExp_MM2(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = p[7]*x[it+nt*7] + p[8];
y[it+nt*1] = (pow(p[7],2))*x[it+nt*3];
    
    return;
}


void dydp_geneExp_MM2(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  memset(dydp,0,sizeof(double)*2*np*nt);
  switch (plist[ip]) {
  case 7: {
dydp[it+nt*(0+ip*ny)] = x[it+nt*7];
dydp[it+nt*(1+ip*ny)] = 2.0*p[7]*x[it+nt*3];

  } break;

  case 8: {
dydp[it+nt*(0+ip*ny)] = 1.0;

  } break;

  }
  }
  
  return;
}


void dydx_geneExp_MM2(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*28);
dydx[7] = (pow(p[7],2));
dydx[14] = p[7];
  
  return;
}


void sy_geneExp_MM2(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)] + x[it+nt*7];
sy[it+nt*(1+np*ny)] = 2.0*p[7]*x[it+nt*3] + (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 8: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)] + 1.0;
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  case 9: {
sy[it+nt*(0+np*ny)] = p[7]*sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = (pow(p[7],2))*sx[it+nt*(3+np*nx)];

  } break;

  }
  
  return;
}
int root_geneExp_MM2(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_geneExp_MM2(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_geneExp_MM2(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_geneExp_MM2(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_geneExp_MM2(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
