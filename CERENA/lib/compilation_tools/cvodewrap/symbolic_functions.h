#include <math.h>

/* ge */

bool cw_ge(double a, double b);
bool Dcw_ge(int id, double a, double b);

/* gt */

bool cw_gt(double a, double b);

bool Dcw_gt(int id, double a, double b);

/* le */

bool cw_le(double a, double b);
bool Dcw_le(int id, double a, double b);

/* lt */

bool cw_lt(double a, double b);

bool Dcw_lt(int id, double a, double b);

/* If */

double cw_if(bool condition, double truepart, double falsepart);
double Dcw_if(int id, bool condition, double truepart, double falsepart);

/* min */

double cw_min(double a, double b);
double Dcw_min(int id,double a, double b);

double cw_max(double a, double b);
double Dcw_max(int id,double a, double b);

/* and */

bool cw_and(bool a, bool b);
double Dcw_and(int id,bool a, bool b);

/* or */

bool cw_or(bool a, bool b);
double Dcw_or(int id,bool a, bool b);


/* step */

double heaviside(double x);
double dirac(double x);
double Ddirac(int ip, double x);

/* splines */

double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);

double Dspline3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double Dspline_pos3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

double Dspline4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double Dspline_pos4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

double Dspline5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double Dspline_pos5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

double Dspline10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double Dspline_pos10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);