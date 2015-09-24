#include "geneExp_LNA.h"
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


 int xdot_geneExp_LNA(realtype t, N_Vector x, N_Vector xdot, void *user_data)
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
xdot_tmp[0] = (p[0]*x_tmp[6] + p[1]*x_tmp[7] - 2.0*k[0]*p[0]*x_tmp[0] + 2.0*k[0]*p[1]*x_tmp[8] - 2.0*(pow(k[0],2))*p[6]*x_tmp[0]*x_tmp[13] - 2.0*(pow(k[0],2))*p[6]*x_tmp[6]*x_tmp[11] + k[0]*p[6]*x_tmp[6]*x_tmp[13])/k[0];
xdot_tmp[1] = (p[0]*x_tmp[6] + p[1]*x_tmp[7] - 2.0*k[0]*p[1]*x_tmp[1] + 2.0*k[0]*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[6]*x_tmp[6]*x_tmp[12] + 2.0*(pow(k[0],2))*p[6]*x_tmp[8]*x_tmp[13] + k[0]*p[6]*x_tmp[6]*x_tmp[13])/k[0];
xdot_tmp[2] = (p[3]*x_tmp[3] + p[2]*x_tmp[7] - 2.0*k[0]*p[3]*x_tmp[2] + 2.0*k[0]*p[2]*x_tmp[10])/k[0];
xdot_tmp[3] = p[2]*x_tmp[7] - 1.0*p[3]*x_tmp[3];
xdot_tmp[4] = (p[4]*x_tmp[3] + p[5]*x_tmp[13] + 2.0*k[0]*p[4]*x_tmp[5] - 2.0*k[0]*p[5]*x_tmp[4])/k[0];
xdot_tmp[5] = p[4]*x_tmp[2] - 1.0*p[3]*x_tmp[5] - 1.0*p[5]*x_tmp[5] + p[2]*x_tmp[12];
xdot_tmp[6] = p[1]*x_tmp[7] - 1.0*p[0]*x_tmp[6] - 1.0*k[0]*p[6]*x_tmp[6]*x_tmp[13];
xdot_tmp[7] = p[0]*x_tmp[6] - 1.0*p[1]*x_tmp[7] + k[0]*p[6]*x_tmp[6]*x_tmp[13];
xdot_tmp[8] = -(1.0*((pow(k[0],2))*p[1]*x_tmp[8] - 1.0*(pow(k[0],2))*p[1]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) + (pow(k[0],2))*x_tmp[8]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[0]*x_tmp[6] + k[0]*p[1]*x_tmp[7] - 1.0*(pow(k[0],3))*p[6]*x_tmp[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[6]*x_tmp[12]))/(pow(k[0],2));
xdot_tmp[9] = p[2]*x_tmp[8] - 1.0*p[0]*x_tmp[9] + p[1]*x_tmp[10] - 1.0*p[3]*x_tmp[9] - 1.0*k[0]*p[6]*x_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*x_tmp[9]*x_tmp[13];
xdot_tmp[10] = p[2]*x_tmp[1] + p[0]*x_tmp[9] - 1.0*p[1]*x_tmp[10] - 1.0*p[3]*x_tmp[10] + k[0]*p[6]*x_tmp[5]*x_tmp[6] + k[0]*p[6]*x_tmp[9]*x_tmp[13];
xdot_tmp[11] = p[1]*x_tmp[12] - 1.0*p[0]*x_tmp[11] + p[4]*x_tmp[9] - 1.0*p[5]*x_tmp[11] - 1.0*k[0]*p[6]*x_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*x_tmp[11]*x_tmp[13];
xdot_tmp[12] = p[0]*x_tmp[11] - 1.0*p[1]*x_tmp[12] + p[4]*x_tmp[10] - 1.0*p[5]*x_tmp[12] + k[0]*p[6]*x_tmp[4]*x_tmp[6] + k[0]*p[6]*x_tmp[11]*x_tmp[13];
xdot_tmp[13] = p[4]*x_tmp[3] - 1.0*p[5]*x_tmp[13];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_geneExp_LNA(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*224);

  for (ixB=0; ixB<224; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_geneExp_LNA(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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


 void x0_geneExp_LNA(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*14);
x0_tmp[3] = -(1.0*(p[9]*(k[3] - 1.0) - 1.0*k[3]*k[17]))/k[0];
x0_tmp[6] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[7] = (k[2]*k[16])/k[0];
x0_tmp[13] = (k[4]*k[18])/k[0];
  
  
  return;
}


 int Jv_geneExp_LNA(N_Vector v, N_Vector Jv, realtype t,
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
Jv_tmp[0] = 2.0*p[1]*v_tmp[8] + (p[1]*v_tmp[7])/k[0] + (v_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*v_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*v_tmp[13])/k[0] - 2.0*k[0]*p[6]*v_tmp[11]*x_tmp[6];
Jv_tmp[1] = (p[1]*v_tmp[7])/k[0] - 2.0*p[1]*v_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*v_tmp[8])/k[0] + (v_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*v_tmp[13])/k[0] + 2.0*k[0]*p[6]*v_tmp[12]*x_tmp[6];
Jv_tmp[2] = 2.0*p[2]*v_tmp[10] - 2.0*p[3]*v_tmp[2] + (p[3]*v_tmp[3])/k[0] + (p[2]*v_tmp[7])/k[0];
Jv_tmp[3] = p[2]*v_tmp[7] - 1.0*p[3]*v_tmp[3];
Jv_tmp[4] = 2.0*p[4]*v_tmp[5] - 2.0*p[5]*v_tmp[4] + (p[4]*v_tmp[3])/k[0] + (p[5]*v_tmp[13])/k[0];
Jv_tmp[5] = p[4]*v_tmp[2] + p[2]*v_tmp[12] - 1.0*v_tmp[5]*(p[3] + p[5]);
Jv_tmp[6] = p[1]*v_tmp[7] - 1.0*v_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*v_tmp[13]*x_tmp[6];
Jv_tmp[7] = v_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*v_tmp[7] + k[0]*p[6]*v_tmp[13]*x_tmp[6];
Jv_tmp[8] = p[1]*v_tmp[1] + v_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*v_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*v_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*v_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*v_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*v_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*v_tmp[12]*x_tmp[6];
Jv_tmp[9] = p[2]*v_tmp[8] - 1.0*v_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*v_tmp[10] - 1.0*k[0]*p[6]*v_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*v_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*v_tmp[13]*x_tmp[9];
Jv_tmp[10] = p[2]*v_tmp[1] - 1.0*v_tmp[10]*(p[1] + p[3]) + v_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*v_tmp[5]*x_tmp[6] + k[0]*p[6]*v_tmp[6]*x_tmp[5] + k[0]*p[6]*v_tmp[13]*x_tmp[9];
Jv_tmp[11] = p[1]*v_tmp[12] - 1.0*v_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*v_tmp[9] - 1.0*k[0]*p[6]*v_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*v_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*v_tmp[13]*x_tmp[11];
Jv_tmp[12] = p[4]*v_tmp[10] - 1.0*v_tmp[12]*(p[1] + p[5]) + v_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*v_tmp[4]*x_tmp[6] + k[0]*p[6]*v_tmp[6]*x_tmp[4] + k[0]*p[6]*v_tmp[13]*x_tmp[11];
Jv_tmp[13] = p[4]*v_tmp[3] - 1.0*p[5]*v_tmp[13];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_geneExp_LNA(N_Vector vB, N_Vector JvB, realtype t,
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
JvB_tmp[0] = ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*vB_tmp[0])/k[0] - 1.0*vB_tmp[8]*(p[0] + k[0]*p[6]*x_tmp[13]);
JvB_tmp[1] = 2.0*p[1]*vB_tmp[1] - 1.0*p[1]*vB_tmp[8] - 1.0*p[2]*vB_tmp[10];
JvB_tmp[2] = 2.0*p[3]*vB_tmp[2] - 1.0*p[4]*vB_tmp[5];
JvB_tmp[3] = p[3]*vB_tmp[3] - 1.0*p[4]*vB_tmp[13] - (1.0*p[3]*vB_tmp[2])/k[0] - (1.0*p[4]*vB_tmp[4])/k[0];
JvB_tmp[4] = 2.0*p[5]*vB_tmp[4] + k[0]*p[6]*vB_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[6];
JvB_tmp[5] = vB_tmp[5]*(p[3] + p[5]) - 2.0*p[4]*vB_tmp[4] + k[0]*p[6]*vB_tmp[9]*x_tmp[6] - 1.0*k[0]*p[6]*vB_tmp[10]*x_tmp[6];
JvB_tmp[6] = vB_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*vB_tmp[7]*(p[0] + k[0]*p[6]*x_tmp[13]) + (vB_tmp[8]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*vB_tmp[0]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*vB_tmp[1]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + k[0]*p[6]*vB_tmp[9]*x_tmp[5] - 1.0*k[0]*p[6]*vB_tmp[10]*x_tmp[5] + k[0]*p[6]*vB_tmp[11]*x_tmp[4] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[4];
JvB_tmp[7] = p[1]*vB_tmp[7] - 1.0*p[1]*vB_tmp[6] - 1.0*p[2]*vB_tmp[3] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[1]*vB_tmp[1])/k[0] - (1.0*p[2]*vB_tmp[2])/k[0] + (p[1]*vB_tmp[8])/k[0];
JvB_tmp[8] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*vB_tmp[8])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[9] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*vB_tmp[1])/k[0] - 2.0*p[1]*vB_tmp[0];
JvB_tmp[9] = vB_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) - 1.0*p[4]*vB_tmp[11] - 1.0*vB_tmp[10]*(p[0] + k[0]*p[6]*x_tmp[13]);
JvB_tmp[10] = vB_tmp[10]*(p[1] + p[3]) - 1.0*p[1]*vB_tmp[9] - 1.0*p[4]*vB_tmp[12] - 2.0*p[2]*vB_tmp[2];
JvB_tmp[11] = vB_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) - 1.0*vB_tmp[12]*(p[0] + k[0]*p[6]*x_tmp[13]) + 2.0*k[0]*p[6]*vB_tmp[0]*x_tmp[6] - 1.0*k[0]*p[6]*vB_tmp[8]*x_tmp[6];
JvB_tmp[12] = vB_tmp[12]*(p[1] + p[5]) - 1.0*p[1]*vB_tmp[11] - 1.0*p[2]*vB_tmp[5] - 2.0*k[0]*p[6]*vB_tmp[1]*x_tmp[6] + k[0]*p[6]*vB_tmp[8]*x_tmp[6];
JvB_tmp[13] = p[5]*vB_tmp[13] + (vB_tmp[8]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[5]*vB_tmp[4])/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*vB_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*vB_tmp[1])/k[0] + k[0]*p[6]*vB_tmp[6]*x_tmp[6] - 1.0*k[0]*p[6]*vB_tmp[7]*x_tmp[6] + k[0]*p[6]*vB_tmp[9]*x_tmp[9] - 1.0*k[0]*p[6]*vB_tmp[10]*x_tmp[9] + k[0]*p[6]*vB_tmp[11]*x_tmp[11] - 1.0*k[0]*p[6]*vB_tmp[12]*x_tmp[11];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_geneExp_LNA(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_geneExp_LNA(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_geneExp_LNA(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
J->data[0] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13]))/k[0];
J->data[8] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[15] = -2.0*p[1];
J->data[22] = p[1];
J->data[24] = p[2];
J->data[30] = -2.0*p[3];
J->data[33] = p[4];
J->data[44] = p[3]/k[0];
J->data[45] = -1.0*p[3];
J->data[46] = p[4]/k[0];
J->data[55] = p[4];
J->data[60] = -2.0*p[5];
J->data[67] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[68] = k[0]*p[6]*x_tmp[6];
J->data[74] = 2.0*p[4];
J->data[75] = - 1.0*p[3] - 1.0*p[5];
J->data[79] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[80] = k[0]*p[6]*x_tmp[6];
J->data[84] = (p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13])/k[0];
J->data[85] = (p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13])/k[0];
J->data[90] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[91] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[92] = -(1.0*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2));
J->data[93] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[94] = k[0]*p[6]*x_tmp[5];
J->data[95] = -1.0*k[0]*p[6]*x_tmp[4];
J->data[96] = k[0]*p[6]*x_tmp[4];
J->data[98] = p[1]/k[0];
J->data[99] = p[1]/k[0];
J->data[100] = p[2]/k[0];
J->data[101] = p[2];
J->data[104] = p[1];
J->data[105] = -1.0*p[1];
J->data[106] = -(1.0*p[1])/k[0];
J->data[112] = 2.0*p[1];
J->data[113] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])/k[0];
J->data[120] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13])))/(pow(k[0],2));
J->data[121] = p[2];
J->data[135] = - 1.0*p[0] - 1.0*p[3] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[136] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[137] = p[4];
J->data[142] = 2.0*p[2];
J->data[149] = p[1];
J->data[150] = - 1.0*p[1] - 1.0*p[3];
J->data[152] = p[4];
J->data[154] = -2.0*k[0]*p[6]*x_tmp[6];
J->data[162] = k[0]*p[6]*x_tmp[6];
J->data[165] = - 1.0*p[0] - 1.0*p[5] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[166] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[169] = 2.0*k[0]*p[6]*x_tmp[6];
J->data[173] = p[2];
J->data[176] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[179] = p[1];
J->data[180] = - 1.0*p[1] - 1.0*p[5];
J->data[182] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6]))/k[0];
J->data[183] = (2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])/k[0];
J->data[186] = p[5]/k[0];
J->data[188] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[189] = k[0]*p[6]*x_tmp[6];
J->data[190] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2));
J->data[191] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[192] = k[0]*p[6]*x_tmp[9];
J->data[193] = -1.0*k[0]*p[6]*x_tmp[11];
J->data[194] = k[0]*p[6]*x_tmp[11];
J->data[195] = -1.0*p[5];

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_geneExp_LNA(realtype t, N_Vector x,
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
  J->rowvals[3] = 8;
  J->rowvals[4] = 10;
  J->rowvals[5] = 2;
  J->rowvals[6] = 5;
  J->rowvals[7] = 2;
  J->rowvals[8] = 3;
  J->rowvals[9] = 4;
  J->rowvals[10] = 13;
  J->rowvals[11] = 4;
  J->rowvals[12] = 11;
  J->rowvals[13] = 12;
  J->rowvals[14] = 4;
  J->rowvals[15] = 5;
  J->rowvals[16] = 9;
  J->rowvals[17] = 10;
  J->rowvals[18] = 0;
  J->rowvals[19] = 1;
  J->rowvals[20] = 6;
  J->rowvals[21] = 7;
  J->rowvals[22] = 8;
  J->rowvals[23] = 9;
  J->rowvals[24] = 10;
  J->rowvals[25] = 11;
  J->rowvals[26] = 12;
  J->rowvals[27] = 0;
  J->rowvals[28] = 1;
  J->rowvals[29] = 2;
  J->rowvals[30] = 3;
  J->rowvals[31] = 6;
  J->rowvals[32] = 7;
  J->rowvals[33] = 8;
  J->rowvals[34] = 0;
  J->rowvals[35] = 1;
  J->rowvals[36] = 8;
  J->rowvals[37] = 9;
  J->rowvals[38] = 9;
  J->rowvals[39] = 10;
  J->rowvals[40] = 11;
  J->rowvals[41] = 2;
  J->rowvals[42] = 9;
  J->rowvals[43] = 10;
  J->rowvals[44] = 12;
  J->rowvals[45] = 0;
  J->rowvals[46] = 8;
  J->rowvals[47] = 11;
  J->rowvals[48] = 12;
  J->rowvals[49] = 1;
  J->rowvals[50] = 5;
  J->rowvals[51] = 8;
  J->rowvals[52] = 11;
  J->rowvals[53] = 12;
  J->rowvals[54] = 0;
  J->rowvals[55] = 1;
  J->rowvals[56] = 4;
  J->rowvals[57] = 6;
  J->rowvals[58] = 7;
  J->rowvals[59] = 8;
  J->rowvals[60] = 9;
  J->rowvals[61] = 10;
  J->rowvals[62] = 11;
  J->rowvals[63] = 12;
  J->rowvals[64] = 13;
  J->colptrs[0] = 0;
  J->colptrs[1] = 2;
  J->colptrs[2] = 5;
  J->colptrs[3] = 7;
  J->colptrs[4] = 11;
  J->colptrs[5] = 14;
  J->colptrs[6] = 18;
  J->colptrs[7] = 27;
  J->colptrs[8] = 34;
  J->colptrs[9] = 38;
  J->colptrs[10] = 41;
  J->colptrs[11] = 45;
  J->colptrs[12] = 49;
  J->colptrs[13] = 54;
  J->colptrs[14] = 65;
J->data[0] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13]))/k[0];
J->data[1] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[2] = -2.0*p[1];
J->data[3] = p[1];
J->data[4] = p[2];
J->data[5] = -2.0*p[3];
J->data[6] = p[4];
J->data[7] = p[3]/k[0];
J->data[8] = -1.0*p[3];
J->data[9] = p[4]/k[0];
J->data[10] = p[4];
J->data[11] = -2.0*p[5];
J->data[12] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[13] = k[0]*p[6]*x_tmp[6];
J->data[14] = 2.0*p[4];
J->data[15] = - 1.0*p[3] - 1.0*p[5];
J->data[16] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[17] = k[0]*p[6]*x_tmp[6];
J->data[18] = (p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13])/k[0];
J->data[19] = (p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13])/k[0];
J->data[20] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[21] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[22] = -(1.0*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2));
J->data[23] = -1.0*k[0]*p[6]*x_tmp[5];
J->data[24] = k[0]*p[6]*x_tmp[5];
J->data[25] = -1.0*k[0]*p[6]*x_tmp[4];
J->data[26] = k[0]*p[6]*x_tmp[4];
J->data[27] = p[1]/k[0];
J->data[28] = p[1]/k[0];
J->data[29] = p[2]/k[0];
J->data[30] = p[2];
J->data[31] = p[1];
J->data[32] = -1.0*p[1];
J->data[33] = -(1.0*p[1])/k[0];
J->data[34] = 2.0*p[1];
J->data[35] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])/k[0];
J->data[36] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13])))/(pow(k[0],2));
J->data[37] = p[2];
J->data[38] = - 1.0*p[0] - 1.0*p[3] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[39] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[40] = p[4];
J->data[41] = 2.0*p[2];
J->data[42] = p[1];
J->data[43] = - 1.0*p[1] - 1.0*p[3];
J->data[44] = p[4];
J->data[45] = -2.0*k[0]*p[6]*x_tmp[6];
J->data[46] = k[0]*p[6]*x_tmp[6];
J->data[47] = - 1.0*p[0] - 1.0*p[5] - 1.0*k[0]*p[6]*x_tmp[13];
J->data[48] = p[0] + k[0]*p[6]*x_tmp[13];
J->data[49] = 2.0*k[0]*p[6]*x_tmp[6];
J->data[50] = p[2];
J->data[51] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[52] = p[1];
J->data[53] = - 1.0*p[1] - 1.0*p[5];
J->data[54] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6]))/k[0];
J->data[55] = (2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])/k[0];
J->data[56] = p[5]/k[0];
J->data[57] = -1.0*k[0]*p[6]*x_tmp[6];
J->data[58] = k[0]*p[6]*x_tmp[6];
J->data[59] = -(1.0*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2));
J->data[60] = -1.0*k[0]*p[6]*x_tmp[9];
J->data[61] = k[0]*p[6]*x_tmp[9];
J->data[62] = -1.0*k[0]*p[6]*x_tmp[11];
J->data[63] = k[0]*p[6]*x_tmp[11];
J->data[64] = -1.0*p[5];
  return(0);
}


 int JBBand_geneExp_LNA(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_geneExp_LNA(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_geneExp_LNA(long int N, realtype t, N_Vector x,
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
JB->data[0] = (2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])/k[0];
JB->data[6] = -(1.0*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0];
JB->data[7] = -(1.0*p[1])/k[0];
JB->data[8] = -2.0*p[1];
JB->data[11] = 2.0*k[0]*p[6]*x_tmp[6];
JB->data[13] = (2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])/k[0];
JB->data[15] = 2.0*p[1];
JB->data[20] = -(1.0*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0];
JB->data[21] = -(1.0*p[1])/k[0];
JB->data[22] = -(1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13]))/k[0];
JB->data[26] = -2.0*k[0]*p[6]*x_tmp[6];
JB->data[27] = -(1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6]))/k[0];
JB->data[30] = 2.0*p[3];
JB->data[31] = -(1.0*p[3])/k[0];
JB->data[35] = -(1.0*p[2])/k[0];
JB->data[38] = -2.0*p[2];
JB->data[45] = p[3];
JB->data[49] = -1.0*p[2];
JB->data[59] = -(1.0*p[4])/k[0];
JB->data[60] = 2.0*p[5];
JB->data[61] = -2.0*p[4];
JB->data[69] = -(1.0*p[5])/k[0];
JB->data[72] = -1.0*p[4];
JB->data[75] = p[3] + p[5];
JB->data[82] = -1.0*p[2];
JB->data[90] = p[0] + k[0]*p[6]*x_tmp[13];
JB->data[91] = -1.0*p[1];
JB->data[97] = k[0]*p[6]*x_tmp[6];
JB->data[104] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
JB->data[105] = p[1];
JB->data[111] = -1.0*k[0]*p[6]*x_tmp[6];
JB->data[112] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
JB->data[113] = -1.0*p[1];
JB->data[118] = (k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12])/(pow(k[0],2));
JB->data[119] = p[1]/k[0];
JB->data[120] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))/(pow(k[0],2));
JB->data[123] = -1.0*k[0]*p[6]*x_tmp[6];
JB->data[124] = k[0]*p[6]*x_tmp[6];
JB->data[125] = ((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8])/(pow(k[0],2));
JB->data[131] = k[0]*p[6]*x_tmp[6];
JB->data[132] = k[0]*p[6]*x_tmp[5];
JB->data[134] = -1.0*p[2];
JB->data[135] = p[0] + p[3] + k[0]*p[6]*x_tmp[13];
JB->data[136] = -1.0*p[1];
JB->data[139] = k[0]*p[6]*x_tmp[9];
JB->data[141] = -1.0*p[2];
JB->data[145] = -1.0*k[0]*p[6]*x_tmp[6];
JB->data[146] = -1.0*k[0]*p[6]*x_tmp[5];
JB->data[149] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
JB->data[150] = p[1] + p[3];
JB->data[153] = -1.0*k[0]*p[6]*x_tmp[9];
JB->data[158] = k[0]*p[6]*x_tmp[6];
JB->data[160] = k[0]*p[6]*x_tmp[4];
JB->data[163] = -1.0*p[4];
JB->data[165] = p[0] + p[5] + k[0]*p[6]*x_tmp[13];
JB->data[166] = -1.0*p[1];
JB->data[167] = k[0]*p[6]*x_tmp[11];
JB->data[172] = -1.0*k[0]*p[6]*x_tmp[6];
JB->data[174] = -1.0*k[0]*p[6]*x_tmp[4];
JB->data[178] = -1.0*p[4];
JB->data[179] = - 1.0*p[0] - 1.0*k[0]*p[6]*x_tmp[13];
JB->data[180] = p[1] + p[5];
JB->data[181] = -1.0*k[0]*p[6]*x_tmp[11];
JB->data[185] = -1.0*p[4];
JB->data[195] = p[5];

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_geneExp_LNA(realtype t, N_Vector x,
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
  JB->rowvals[1] = 6;
  JB->rowvals[2] = 7;
  JB->rowvals[3] = 8;
  JB->rowvals[4] = 11;
  JB->rowvals[5] = 13;
  JB->rowvals[6] = 1;
  JB->rowvals[7] = 6;
  JB->rowvals[8] = 7;
  JB->rowvals[9] = 8;
  JB->rowvals[10] = 12;
  JB->rowvals[11] = 13;
  JB->rowvals[12] = 2;
  JB->rowvals[13] = 3;
  JB->rowvals[14] = 7;
  JB->rowvals[15] = 10;
  JB->rowvals[16] = 3;
  JB->rowvals[17] = 7;
  JB->rowvals[18] = 3;
  JB->rowvals[19] = 4;
  JB->rowvals[20] = 5;
  JB->rowvals[21] = 13;
  JB->rowvals[22] = 2;
  JB->rowvals[23] = 5;
  JB->rowvals[24] = 12;
  JB->rowvals[25] = 6;
  JB->rowvals[26] = 7;
  JB->rowvals[27] = 13;
  JB->rowvals[28] = 6;
  JB->rowvals[29] = 7;
  JB->rowvals[30] = 13;
  JB->rowvals[31] = 0;
  JB->rowvals[32] = 1;
  JB->rowvals[33] = 6;
  JB->rowvals[34] = 7;
  JB->rowvals[35] = 8;
  JB->rowvals[36] = 11;
  JB->rowvals[37] = 12;
  JB->rowvals[38] = 13;
  JB->rowvals[39] = 5;
  JB->rowvals[40] = 6;
  JB->rowvals[41] = 8;
  JB->rowvals[42] = 9;
  JB->rowvals[43] = 10;
  JB->rowvals[44] = 13;
  JB->rowvals[45] = 1;
  JB->rowvals[46] = 5;
  JB->rowvals[47] = 6;
  JB->rowvals[48] = 9;
  JB->rowvals[49] = 10;
  JB->rowvals[50] = 13;
  JB->rowvals[51] = 4;
  JB->rowvals[52] = 6;
  JB->rowvals[53] = 9;
  JB->rowvals[54] = 11;
  JB->rowvals[55] = 12;
  JB->rowvals[56] = 13;
  JB->rowvals[57] = 4;
  JB->rowvals[58] = 6;
  JB->rowvals[59] = 10;
  JB->rowvals[60] = 11;
  JB->rowvals[61] = 12;
  JB->rowvals[62] = 13;
  JB->rowvals[63] = 3;
  JB->rowvals[64] = 13;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 6;
  JB->colptrs[2] = 12;
  JB->colptrs[3] = 16;
  JB->colptrs[4] = 18;
  JB->colptrs[5] = 22;
  JB->colptrs[6] = 25;
  JB->colptrs[7] = 28;
  JB->colptrs[8] = 31;
  JB->colptrs[9] = 39;
  JB->colptrs[10] = 45;
  JB->colptrs[11] = 51;
  JB->colptrs[12] = 57;
  JB->colptrs[13] = 63;
  JB->colptrs[14] = 65;
  return(0);
}


 int sx_geneExp_LNA(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] - (1.0*(2.0*k[0]*x_tmp[0] - 1.0*x_tmp[6]))/k[0] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (2.0*k[0]*x_tmp[8] + x_tmp[6])/k[0] - 2.0*p[1]*sx_tmp[1] + (p[1]*sx_tmp[7])/k[0] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + x_tmp[6] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + x_tmp[9] + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*x_tmp[11] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + x_tmp[11] + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 1: {
sxdot_tmp[0] = (2.0*k[0]*x_tmp[8] + x_tmp[7])/k[0] + 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - (1.0*(2.0*k[0]*x_tmp[1] - 1.0*x_tmp[7]))/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) + x_tmp[7] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] - 1.0*x_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] + x_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*x_tmp[10] + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] + x_tmp[12] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*x_tmp[12] + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 2: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = (2.0*k[0]*x_tmp[10] + x_tmp[7])/k[0] - 2.0*p[3]*sx_tmp[2] + 2.0*p[2]*sx_tmp[10] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3] + x_tmp[7];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]) + x_tmp[12];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] + x_tmp[8] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + x_tmp[1] + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 3: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] - (1.0*(2.0*k[0]*x_tmp[2] - 1.0*x_tmp[3]))/k[0] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3] - 1.0*x_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]) - 1.0*x_tmp[5];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*x_tmp[10] + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 4: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = (2.0*k[0]*x_tmp[5] + x_tmp[3])/k[0] + 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]) + x_tmp[2];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] + x_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + x_tmp[10] + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13] + x_tmp[3];

  } break;

  case 5: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] - (1.0*(2.0*k[0]*x_tmp[4] - 1.0*x_tmp[13]))/k[0] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]) - 1.0*x_tmp[5];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*x_tmp[11] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*x_tmp[12] + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13] - 1.0*x_tmp[13];

  } break;

  case 6: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] - (1.0*(2.0*(pow(k[0],2))*x_tmp[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[6]*x_tmp[11] - 1.0*k[0]*x_tmp[6]*x_tmp[13]))/k[0] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (2.0*(pow(k[0],2))*x_tmp[6]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[8]*x_tmp[13] + k[0]*x_tmp[6]*x_tmp[13])/k[0] - 2.0*p[1]*sx_tmp[1] + (p[1]*sx_tmp[7])/k[0] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*x_tmp[6]*x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*x_tmp[6]*x_tmp[13] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*((pow(k[0],2))*x_tmp[6]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[0]*x_tmp[13] + (pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[8]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*x_tmp[5]*x_tmp[6] - 1.0*k[0]*x_tmp[9]*x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*x_tmp[5]*x_tmp[6] + k[0]*x_tmp[9]*x_tmp[13] + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*x_tmp[4]*x_tmp[6] - 1.0*k[0]*x_tmp[11]*x_tmp[13] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*x_tmp[4]*x_tmp[6] + k[0]*x_tmp[11]*x_tmp[13] + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 7: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 8: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  case 9: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[8] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[6]*(p[0] - 2.0*(pow(k[0],2))*p[6]*x_tmp[11] + k[0]*p[6]*x_tmp[13]))/k[0] - (1.0*(2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[0])/k[0] - (1.0*(2.0*(pow(k[0],2))*p[6]*x_tmp[0] - 1.0*k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] - 2.0*k[0]*p[6]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[1] = (p[1]*sx_tmp[7])/k[0] - 2.0*p[1]*sx_tmp[1] + ((2.0*k[0]*p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[13])*sx_tmp[8])/k[0] + (sx_tmp[6]*(p[0] + 2.0*(pow(k[0],2))*p[6]*x_tmp[12] + k[0]*p[6]*x_tmp[13]))/k[0] + ((2.0*(pow(k[0],2))*p[6]*x_tmp[8] + k[0]*p[6]*x_tmp[6])*sx_tmp[13])/k[0] + 2.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[10] - 2.0*p[3]*sx_tmp[2] + (p[3]*sx_tmp[3])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[3] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[3];
sxdot_tmp[4] = 2.0*p[4]*sx_tmp[5] - 2.0*p[5]*sx_tmp[4] + (p[4]*sx_tmp[3])/k[0] + (p[5]*sx_tmp[13])/k[0];
sxdot_tmp[5] = p[4]*sx_tmp[2] + p[2]*sx_tmp[12] - 1.0*sx_tmp[5]*(p[3] + p[5]);
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[7] = sx_tmp[6]*(p[0] + k[0]*p[6]*x_tmp[13]) - 1.0*p[1]*sx_tmp[7] + k[0]*p[6]*sx_tmp[13]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[1] + sx_tmp[0]*(p[0] + k[0]*p[6]*x_tmp[13]) - (1.0*sx_tmp[6]*(k[0]*p[0] - 1.0*(pow(k[0],3))*p[6]*x_tmp[11] + (pow(k[0],2))*p[6]*x_tmp[13] + (pow(k[0],3))*p[6]*x_tmp[12]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*((pow(k[0],2))*p[6]*x_tmp[6] - 1.0*(pow(k[0],3))*p[6]*x_tmp[0] + (pow(k[0],3))*p[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*(p[0] + k[0]*p[6]*x_tmp[13]))*sx_tmp[8])/(pow(k[0],2)) + k[0]*p[6]*sx_tmp[11]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[9] = p[2]*sx_tmp[8] - 1.0*sx_tmp[9]*(p[0] + p[3] + k[0]*p[6]*x_tmp[13]) + p[1]*sx_tmp[10] - 1.0*k[0]*p[6]*sx_tmp[5]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[5] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[10] = p[2]*sx_tmp[1] - 1.0*sx_tmp[10]*(p[1] + p[3]) + sx_tmp[9]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[5]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[5] + k[0]*p[6]*sx_tmp[13]*x_tmp[9];
sxdot_tmp[11] = p[1]*sx_tmp[12] - 1.0*sx_tmp[11]*(p[0] + p[5] + k[0]*p[6]*x_tmp[13]) + p[4]*sx_tmp[9] - 1.0*k[0]*p[6]*sx_tmp[4]*x_tmp[6] - 1.0*k[0]*p[6]*sx_tmp[6]*x_tmp[4] - 1.0*k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[12] = p[4]*sx_tmp[10] - 1.0*sx_tmp[12]*(p[1] + p[5]) + sx_tmp[11]*(p[0] + k[0]*p[6]*x_tmp[13]) + k[0]*p[6]*sx_tmp[4]*x_tmp[6] + k[0]*p[6]*sx_tmp[6]*x_tmp[4] + k[0]*p[6]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = p[4]*sx_tmp[3] - 1.0*p[5]*sx_tmp[13];

  } break;

  }
 for (ix=0; ix<14; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_geneExp_LNA(int ip, N_Vector sx0, void *user_data)
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
sx0_tmp[3] = -(1.0*(k[3] - 1.0))/k[0];

  } break;

  }

  return;
}


void y_geneExp_LNA(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*6];
y[it+nt*1] = x[it+nt*7];
y[it+nt*2] = x[it+nt*3];
y[it+nt*3] = x[it+nt*13];
y[it+nt*4] = x[it+nt*0];
y[it+nt*5] = x[it+nt*8];
y[it+nt*6] = x[it+nt*1];
y[it+nt*7] = x[it+nt*9];
y[it+nt*8] = x[it+nt*10];
y[it+nt*9] = x[it+nt*2];
y[it+nt*10] = x[it+nt*11];
y[it+nt*11] = x[it+nt*12];
y[it+nt*12] = x[it+nt*5];
y[it+nt*13] = x[it+nt*4];
y[it+nt*14] = p[7]*x[it+nt*13] + p[8];
y[it+nt*15] = (pow(p[7],2))*x[it+nt*4];
    
    return;
}


void dydp_geneExp_LNA(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  memset(dydp,0,sizeof(double)*16*np*nt);
  switch (plist[ip]) {
  case 7: {
dydp[it+nt*(14+ip*ny)] = x[it+nt*13];
dydp[it+nt*(15+ip*ny)] = 2.0*p[7]*x[it+nt*4];

  } break;

  case 8: {
dydp[it+nt*(14+ip*ny)] = 1.0;

  } break;

  }
  }
  
  return;
}


void dydx_geneExp_LNA(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*224);
dydx[4] = 1.0;
dydx[22] = 1.0;
dydx[41] = 1.0;
dydx[50] = 1.0;
dydx[77] = 1.0;
dydx[79] = (pow(p[7],2));
dydx[92] = 1.0;
dydx[96] = 1.0;
dydx[113] = 1.0;
dydx[133] = 1.0;
dydx[151] = 1.0;
dydx[168] = 1.0;
dydx[186] = 1.0;
dydx[203] = 1.0;
dydx[211] = 1.0;
dydx[222] = p[7];
  
  return;
}


void sy_geneExp_LNA(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 6: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 7: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)] + x[it+nt*13];
sy[it+nt*(15+np*ny)] = 2.0*p[7]*x[it+nt*4] + (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 8: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)] + 1.0;
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  case 9: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(2+np*ny)] = sx[it+nt*(3+np*nx)];
sy[it+nt*(3+np*ny)] = sx[it+nt*(13+np*nx)];
sy[it+nt*(4+np*ny)] = sx[it+nt*(0+np*nx)];
sy[it+nt*(5+np*ny)] = sx[it+nt*(8+np*nx)];
sy[it+nt*(6+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(7+np*ny)] = sx[it+nt*(9+np*nx)];
sy[it+nt*(8+np*ny)] = sx[it+nt*(10+np*nx)];
sy[it+nt*(9+np*ny)] = sx[it+nt*(2+np*nx)];
sy[it+nt*(10+np*ny)] = sx[it+nt*(11+np*nx)];
sy[it+nt*(11+np*ny)] = sx[it+nt*(12+np*nx)];
sy[it+nt*(12+np*ny)] = sx[it+nt*(5+np*nx)];
sy[it+nt*(13+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(14+np*ny)] = p[7]*sx[it+nt*(13+np*nx)];
sy[it+nt*(15+np*ny)] = (pow(p[7],2))*sx[it+nt*(4+np*nx)];

  } break;

  }
  
  return;
}
int root_geneExp_LNA(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_geneExp_LNA(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_geneExp_LNA(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_geneExp_LNA(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_geneExp_LNA(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
