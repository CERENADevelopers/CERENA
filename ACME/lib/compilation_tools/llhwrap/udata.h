#ifndef _MY_UDATA
#define _MY_UDATA

typedef struct {
	double *qpositivex;
	double *u;
    int    *plist;
    int    *np;
    int    *ny;
	double *p;
    double *k;
	} *UserData;
#endif /* _MY_UDATA */
