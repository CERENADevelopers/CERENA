#ifndef _MY_UDATA
#define _MY_UDATA

typedef struct {
	double *qpositivex;
	double *u;
	double *su;
    int    *plist;
    int    *nps;
    int    *nobs;
    double *par;
	double *p;
    double *k;
	double *v;
	double *sv;
	double *dvdx;
    double *dvddx;
	double *dvdu;
	double *dvdp;
	} *UserData;
#endif /* _MY_UDATA */
