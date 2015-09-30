/*
 *  MATLAB usage: [status,t,x,y,sx,sy] = simulate_"model"(tout,theta,kappa,options)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <cvodes/cvodes.h>

#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_klu.h>
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_sparse.h> /* def. of type sparse stuff */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_config.h>  /* ? */
/* #if SUNDIALS_BLAS_LAPACK
 #include <cvodes/cvodes_lapack.h>
 #include <sundials/sundials_lapack.h>
 #endif */
#include <klu.h>
#include <amd.h>
#include <colamd.h>
#include <btf.h>
/*#include <Accelerate/Accelerate.h>*/


/* user functions */
#include "cvodewrapfunctions.c"

/* Accessor macros */
#define Ith(v, i)     NV_Ith_S(v, i-1)        /* i-th vector component i=1..nx */

#define LW_DENSE       1
#define LW_BAND        2
#define LW_LAPACKDENSE 3
#define LW_LAPACKBAND  4
#define LW_DIAG        5
#define LW_SPGMR       6
#define LW_SPBCG       7
#define LW_SPTFQMR     8
#define LW_KLU         9

#define LW_DENSE       1
#define LW_BAND        2

#define LW_NORMAL      1
#define LW_LOGNORMAL   2
#define LW_ONEOUTPUT   5

#define pi             3.141592653589793


int sensi;
double cvodes_atol,cvodes_rtol;
int cvodes_maxsteps;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int cv_status;
    double *status;
    int np, nx, ny, nt, nr, nnz, ndisc; /* integers with problem dimensionality */
    int nmyt, nmyy, nysigmat, nysigmay; /* integers with problem dimensionality */
    int nmtt, nmty, ntsigmat, ntsigmay; /* integers with problem dimensionality */
    int ip, ix, jx, iy, it, ir, idisc; /* integers for indexing in loops */
    int jp; /* integers for indexing in loops */
    int ncheck; /* the number of (internal) checkpoints stored so far */
    int nroots; /* number of found roots */
    
    void *cvode_mem, *cvode_memB; /* pointer to cvodes memory block */
    UserData data; /* user data */
    
    realtype t; /* current time */
    double tstart; /* starting time */
    N_Vector x; /* state vector */
    N_Vector xdot; /* differential state vector */
    N_Vector xB; /* state vector for backwards problem */
    N_Vector xQB; /* state vector for backwards integration problem */
    N_Vector atols_ss;
    N_Vector *sx; /* sensitivity state vector array */
    
    double *ts, *tdata; /* timepoints */
    double *xdotdata; /* returned vector containing differential states*/
    double *Jdata;
    double *rootdata;
    int    *rootidx;
    double *rootvaldata;
    double *rootvaltmp; /* returned vector containing states and sensitivities */
    double *xdata; /* returned vector containing states and sensitivities */
    double *ydata; /* returned vector containing outputs and sensitivities */
    double *llhdata, *chi2data, *llhSdata, *llhS2data;
    double g, *dgdp, *dgdx;
    double r, *drdp, *drdx;
    double rval, *drvaldp, *drvaldx;
    double *my, *ysigma;
    double *mt, *tsigma;
    
    int nmaxroot, nmaxdisc; /* maximum number of roots to record */
    int *rootsfound; /* array of length nr with the indices of the user functions gi found to have a root. For i = 0, . . . ,nr?1, rootsfound[i]?= 0 if gi has a root, and = 0 if not.
                      */
    int ubw, lbw; /* upper and lower bandwith of the jacobian */
    
    long int numsteps;
    long int numrhsevals;
    long int numlinsolvsetups;
    long int numerrtestfails;
    int order;
    long int numnonlinsolviters;
    long int numjacevals;
    long int numliniters;
    long int numconvfails;
    long int numprecevals;
    long int numprecsolves;
    long int numjtimesevals;
    
    double *numstepsdata, *numstepsSdata;
    double *numrhsevalsdata, *numrhsevalsSdata;
    double *numlinsolvsetupsdata, *numlinsolvsetupsSdata;
    double *numerrtestfailsdata, *numerrtestfailsSdata;
    double *orderdata, *orderSdata;
    double *numnonlinsolvitersdata, *numnonlinsolvitersSdata;
    double *numjacevalsdata, *numjacevalsSdata;
    double *numlinitersdata, *numlinitersSdata;
    double *numconvfailsdata, *numconvfailsSdata;
    double *numprecevalsdata, *numprecevalsSdata;
    double *numprecsolvesdata, *numprecsolvesSdata;
    double *numjtimesevalsdata, *numjtimesevalsSdata;
    
    
    double *llhS0;
    realtype *x_tmp, *sx_tmp; /* temporary storage for states and sensitivities */
    realtype *xdot_tmp; /* temporary storage for states and sensitivities */
    realtype *xB_tmp, *xQB_tmp; /* temporary storage for backward states */
    DlsMat Jtmp;
    double *dydx, *dydp, *xnow, *ynow; /* temporary storage for initiation of backward states */
    double *dtdx, *dtdp; /* temporary storage for initiation of backward states */
    double *p;  /* parameter vector */
    int    *plist;  /* index set of parameters for which to compute sensitivities */
    double *plistdata; /* input for plist */
    double *pbar;
    double *xbar;
    
    long int    Nd;     /* number of integration steps between two consecutive checkpoints. */
    
    int ism; /* a flag used to select the sensitivity solution method. Its value can be CV SIMULTANEOUS or CV STAGGERED */
    int linsol; /* linear solver specification */
    int interpType; /* CV_POLYNOMIAL or CV_HERMITE */
    int lmm; /* FORWARD PROBLEM: specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.  */
    int iter; /* FORWARD PROBLEM: specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.  */
    int lmmB; /* BACKWARD PROBLEM: specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.  */
    int iterB; /* BACKWARD PROBLEM: specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.  */
    int which = 0; /* BACKWARD PROBLEM: contains the identifier assigned by cvodes for the newly created back- ward problem. Any call to CVode*B functions requires such an identifier. */
    bool stldet; /* flag controlling stability limit detection */
    double stldetdata;
    bool error_corr = TRUE;
    int pretype= PREC_NONE; /* specifies the type of preconditioning and must be one of: PREC NONE, PREC LEFT, PREC RIGHT, or PREC BOTH. */
    int data_model, event_model;
    int ordering; /* specification of ordering */
    
    double *discs; /* vector containing the time-points of discontinuities */
    double *irdiscs; /* vector containing the index of discontinuities */
    bool rootflag, discflag;
    
    
    
    
    /*
     * ----------------------------------------------------------------
     * EXTRACT INPUTS
     * ----------------------------------------------------------------
     */
    
    /* User data structure */
    data = (UserData) mxMalloc(sizeof *data);
    if(data == NULL) goto error_return;
    
    
    /* solution struct */
    
    if (!prhs[0]) {
        mexErrMsgTxt("No solution struct provided!");
    }
    
    if(mxGetField(prhs[0], 0 ,"status")) { status = mxGetPr(mxGetField(prhs[0], 0 ,"status")); } else { mexErrMsgTxt("Parameter status not specified as field in solution struct!"); }
    
    
    /* time vector */
    
    if (!prhs[1]) {
        mexErrMsgTxt("No time vector provided!");
    }
    
    ts = mxGetPr(prhs[1]);
    
    nt = mxGetM(prhs[1]) * mxGetN(prhs[1]);
    
    /* parameters */
    
    if (!prhs[2]) {
        mexErrMsgTxt("No parameter vector provided!");
    }
    data->p = mxGetPr(prhs[2]);
    p = data->p;
    
    /* constants */
    
    if (!prhs[3]) {
        mexErrMsgTxt("No constant vector provided!");
    }
    data->k = mxGetPr(prhs[3]);
    
    /* data */
    
    if(mxGetField(prhs[5], 0 ,"data_model")) { data_model = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"data_model")); } else { mexErrMsgTxt("Parameter data_model not specified as field in options struct!"); }
    
    /* options */
    if (!prhs[5]) {
        mexErrMsgTxt("No options provided!");
    }
    
    if(mxGetField(prhs[5], 0 ,"tstart")) { tstart = mxGetScalar(mxGetField(prhs[5], 0 ,"tstart")); } else { mexErrMsgTxt("Parameter tstart not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"cvodes_atol")) { cvodes_atol = mxGetScalar(mxGetField(prhs[5], 0 ,"cvodes_atol")); } else { mexErrMsgTxt("Parameter cvodes_atol not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"cvodes_rtol")) { cvodes_rtol = mxGetScalar(mxGetField(prhs[5], 0 ,"cvodes_rtol")); } else { mexErrMsgTxt("Parameter cvodes_rtol not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"cvodes_maxsteps")) { cvodes_maxsteps = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"cvodes_maxsteps")); } else { mexErrMsgTxt("Parameter cvodes_maxsteps not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"nx")) { nx = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"nx")); } else { mexErrMsgTxt("Parameter nx not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"ny")) { ny = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"ny")); } else { mexErrMsgTxt("Parameter ny not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"np")) { np = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"np")); } else { mexErrMsgTxt("Parameter np not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"nr")) { nr = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"nr")); } else { mexErrMsgTxt("Parameter nr not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"ndisc")) { ndisc = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"ndisc")); } else { mexErrMsgTxt("Parameter ndisc not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"nnz")) { nnz = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"nnz")); } else { mexErrMsgTxt("Parameter nnz not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"lmm")) { lmm = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"lmm")); } else {  mexErrMsgTxt("Parameter lmm not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"iter")) { iter = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"iter")); } else { mexErrMsgTxt("Parameter iter not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"lmmB")) { lmmB = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"lmmB")); } else { mexErrMsgTxt("Parameter lmmB not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"iterB")) { iterB = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"iterB")); } else { mexErrMsgTxt("Parameter iterB not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"Nd"))  { Nd = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"Nd")); } else { mexErrMsgTxt("Parameter Nd not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"interpType"))  { interpType = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"interpType")); } else { mexErrMsgTxt("Parameter interpType not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"linsol")) { linsol = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"linsol")); } else { mexErrMsgTxt("Parameter linsol not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"stldet")) { stldetdata = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"stldet")); } else { mexErrMsgTxt("Parameter stldetdata not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"sensi")) { sensi = (int) mxGetScalar(mxGetField(prhs[5], 0 ,"sensi")); } else { mexErrMsgTxt("Parameter sensi not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"ism")) { ism = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"ism")); } else { mexErrMsgTxt("Parameter ism not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"nmaxroot")) { nmaxroot = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"nmaxroot")); } else { mexErrMsgTxt("Parameter nmaxroot not specified as field in options struct!"); }
    if(mxGetField(prhs[5], 0 ,"nmaxdisc")) { nmaxdisc = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"nmaxdisc")); } else { mexErrMsgTxt("Parameter nmaxdisc not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"ubw")) { ubw = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"ubw")); } else { mexErrMsgTxt("Parameter ubw not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"lbw")) { lbw = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"lbw")); } else { mexErrMsgTxt("Parameter lbw not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"event_model")) { event_model = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"event_model")); } else { mexErrMsgTxt("Parameter event_model not specified as field in options struct!"); }
    
    if(mxGetField(prhs[5], 0 ,"ordering")) { ordering = (int)mxGetScalar(mxGetField(prhs[5], 0 ,"ordering")); } else { mexErrMsgTxt("Parameter ordering not specified as field in options struct!"); }
    
    
    if ((int)stldetdata>0.5) {
        stldet = TRUE;
    } else {
        stldet = FALSE;
    }
    
    
    
    if (data_model == LW_ONEOUTPUT) {
        if ( ny>1 | nt>1 ) {
            mexErrMsgTxt("Data model LW_ONEOUTPUT not allowed for more than one time-point or more than one observable!");
        }
    } else {
        
        if (!prhs[4]) {
            mexErrMsgTxt("No data provided!");
        }
        if (mxGetField(prhs[4], 0 ,"Y")) {
            my = mxGetPr(mxGetField(prhs[4], 0 ,"Y"));
            nmyy = mxGetN(mxGetField(prhs[4], 0 ,"Y"));
            nmyt = mxGetM(mxGetField(prhs[4], 0 ,"Y"));
        } else {
            mexErrMsgTxt("Field Y not specified as field in data struct!");
        }
        
        if (mxGetField(prhs[4], 0 ,"Sigma_Y")) {
            ysigma = mxGetPr(mxGetField(prhs[4], 0 ,"Sigma_Y"));
            nysigmay = mxGetN(mxGetField(prhs[4], 0 ,"Sigma_Y"));
            nysigmat = mxGetM(mxGetField(prhs[4], 0 ,"Sigma_Y"));
        } else {
            mexErrMsgTxt("Field Sigma_Y not specified as field in data struct!");
        }
        if (mxGetField(prhs[4], 0 ,"T")) {
            mt = mxGetPr(mxGetField(prhs[4], 0 ,"T"));
            nmty = mxGetN(mxGetField(prhs[4], 0 ,"T"));
            nmtt = mxGetM(mxGetField(prhs[4], 0 ,"T"));
        } else {
            mexErrMsgTxt("Field T not specified as field in data struct!");
        }
        
        if (mxGetField(prhs[4], 0 ,"Sigma_T")) {
            tsigma = mxGetPr(mxGetField(prhs[4], 0 ,"Sigma_T"));
            ntsigmay = mxGetN(mxGetField(prhs[4], 0 ,"Sigma_T"));
            ntsigmat = mxGetM(mxGetField(prhs[4], 0 ,"Sigma_T"));
        } else {
            mexErrMsgTxt("Field Sigma_T not specified as field in data struct!");
        }
        
        if (nmyt != nt) {
            mexErrMsgTxt("Number of time-points in data matrix does not match provided time vector");
        }
        
        if (nysigmat != nt) {
            mexErrMsgTxt("Number of time-points in data-sigma matrix does not match provided time vector");
        }
        
        if (nmyy != ny) {
            mexErrMsgTxt("Number of observables in data matrix does not match provided ny");
        }
        
        if (nysigmay != ny) {
            mexErrMsgTxt("Number of observables in data-sigma matrix does not match provided ny");
        }
        
        if (nmtt != nmaxroot) {
            mexErrMsgTxt("Number of time-points in event matrix does not match provided nmaxroot");
        }
        
        if (ntsigmat != nmaxroot) {
            mexErrMsgTxt("Number of time-points in event-sigma matrix does not match provided nmaxroot");
        }
        
        if (nmty != nr) {
            mexErrMsgTxt("Number of events in event matrix does not match provided nr");
        }
        
        if (ntsigmay != nr) {
            mexErrMsgTxt("Number of events in event-sigma matrix does not match provided nr");
        }
        
    }
    
    
    
    data->np = &RCONST(np);
    data->ny = &RCONST(ny);
    if(mxGetField(prhs[5], 0 ,"qPositiveX")) { data->qpositivex = mxGetData(mxGetField(prhs[5], 0, "qPositiveX")); } else { mexErrMsgTxt("Parameter qPositiveX not specified as field in options struct!"); }
    
    /* plist */
    if (!prhs[6]) {
        mexErrMsgTxt("No parameter list provided!");
    }
    
    if(prhs[6]) {
        plistdata = mxGetPr(prhs[6]);
    }
    
    plist = mxMalloc(np*sizeof(int));
    for (ip=0; ip<np; ip++) {
        plist[ip] = (int)plistdata[ip];
    }
    
    data->plist = plist;
    
    /* pbar */
    if (!prhs[7]) {
        mexErrMsgTxt("No parameter scales provided!");
    }
    
    pbar = mxGetPr(prhs[7]);
    
    /* xscale */
    if (!prhs[8]) {
        mexErrMsgTxt("No state scales provided!");
    }
    
    xbar = mxGetPr(prhs[8]);
    
    
    
    
    /*
     * ----------------------------------------------------------------
     * CVODEINIT
     * ----------------------------------------------------------------
     */
    
    
    
    if (nx > 0) {
        
        /* initialise log-likelihood */
        g = 0; /* data part */
        r = 0; /* event part */
        rval = 0;
        
        /* write initial conditions */
        x = N_VNew_Serial(nx);
        xdot = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);
        
        if(nr>0) rootvaltmp = mxMalloc(nr*sizeof(double));
        if(nr+ndisc>0) rootsfound = mxMalloc((nr+ndisc)*sizeof(double));
        if(nr>0) rootidx = mxMalloc(nr*sizeof(int));
        if(ndisc>0) discs = mxMalloc(nmaxdisc*sizeof(double));
        if(ndisc>0) irdiscs = mxMalloc(nmaxdisc*sizeof(double));
        
        if (x == NULL) goto error_return;
        fx0(x, data);
        
        /* Create CVODES object */
        if (lmm>2||lmm<1) {
            mexErrMsgTxt("Illegal value for lmm!");
        }
        if (iter>2||iter<1) {
            mexErrMsgTxt("Illegal value for iter!");
        }
        cvode_mem = CVodeCreate(lmm, iter);
        if (cvode_mem == NULL) goto error_return;
        
        /* Initialize CVODES solver*/
        *status = cvodewrap_init(cvode_mem, x, tstart);
        if (*status != CV_SUCCESS) goto error_return;
        
        /* Specify integration tolerances */
        *status = CVodeSStolerances(cvode_mem, RCONST(cvodes_rtol), RCONST(cvodes_atol));
        if(*status != CV_SUCCESS) goto error_return;
        
        /* Set optional inputs */
        *status = CVodeSetErrHandlerFn(cvode_mem, cvodewrap_ErrHandlerFn, NULL);
        if(*status != CV_SUCCESS) goto error_return;
        
        *status = CVodeSetUserData(cvode_mem, data); /* attaches userdata*/
        if(*status != CV_SUCCESS) goto error_return;
        
        *status = CVodeSetMaxNumSteps(cvode_mem, cvodes_maxsteps); /* specify maximal number of steps */
        if(*status != CV_SUCCESS) goto error_return;
        
        *status = CVodeSetStabLimDet(cvode_mem, stldet); /* activates stability limit detection */
        if(*status != CV_SUCCESS) goto error_return;
        
        if (nr+ndisc>0) {
            *status = cvodewrap_RootInit(cvode_mem, nr+ndisc); /* activates root detection */
            if(*status != CV_SUCCESS) goto error_return;
        }
        
        /* Attach linear solver module */
        switch (linsol) {
                
                /* DIRECT SOLVERS */
                
            case LW_DENSE:
                *status = CVDense(cvode_mem, nx);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetDenseJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case LW_BAND:
                *status = CVBand(cvode_mem, nx, ubw, lbw);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetBandJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case LW_LAPACKDENSE:
                
                /*#if SUNDIALS_BLAS_LAPACK
                 *status = CVLapackDense(cvode_mem, nx);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 *status = cvodewrap_SetDenseJacFn(cvode_mem);
                 if (*status != CV_SUCCESS) goto error_return;
                 #else*/
                mexErrMsgTxt("Solver currently not supported!");
                /*#endif*/
                break;
                
            case LW_LAPACKBAND:
                
                
                /*#if SUNDIALS_BLAS_LAPACK
                 *status = CVLapackBand(cvode_mem, nx, ubw, lbw);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 *status = cvodewrap_SetBandJacFn(cvode_mem);
                 if (*status != CV_SUCCESS) goto error_return;
                 #else*/
                mexErrMsgTxt("Solver currently not supported!");
                /*#endif*/
                break;
                
            case LW_DIAG:
                *status = CVDiag(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
                /* ITERATIVE SOLVERS */
                
            case LW_SPGMR:
                *status = CVSpgmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case LW_SPBCG:
                *status = CVSpbcg(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case LW_SPTFQMR:
                *status = CVSptfqmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case LW_KLU:
                *status = CVKLU(cvode_mem, nx, nnz);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetSparseJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVKLUSetOrdering(cvode_mem, ordering);
                if (*status != CV_SUCCESS) goto error_return;
                
                
                break;
                
            default:
                mexErrMsgTxt("Invalid choice of solver!");
                break;
        }
        
        
        
    }
    
    /*
     * ----------------------------------------------------------------
     * CVODESENSINIT
     * ----------------------------------------------------------------
     */
    
    if ( sensi >= 1) {
        if(nx>0) {
            /* Allocate space for the adjoint computation */
            *status = CVodeAdjInit(cvode_mem, Nd, interpType);
            if (*status != CV_SUCCESS) goto error_return;
            
            xnow = mxMalloc(nx*sizeof(double));
            ynow = mxMalloc(ny*sizeof(double));
            dydx = mxMalloc(ny*nx*sizeof(double));
            memset(dydx,0,sizeof(double)*ny*nx);
            dydp = mxMalloc(ny*np*sizeof(double));
            memset(dydp,0,sizeof(double)*ny*np);
            llhS0 = mxMalloc(np*sizeof(double));
            dgdp = mxMalloc(np*sizeof(double));
            memset(dgdp,0,sizeof(double)*np);
            dgdx = mxMalloc(nx*nt*sizeof(double));
            memset(dgdx,0,sizeof(double)*nx*nt);
            dtdp = mxMalloc(nr*np*sizeof(double));
            memset(dgdp,0,sizeof(double)*nr*np);
            dtdx = mxMalloc(nr*nx*sizeof(double));
            memset(dgdx,0,sizeof(double)*nr*nx);
            drdp = mxMalloc(np*sizeof(double));
            memset(drdp,0,sizeof(double)*np);
            drdx = mxMalloc(nx*sizeof(double));
            memset(drdx,0,sizeof(double)*nx*nmaxroot);
            drvaldp = mxMalloc(np*sizeof(double));
            memset(drvaldp,0,sizeof(double)*np);
            drvaldx = mxMalloc(nx*nmaxroot*sizeof(double));
            memset(drvaldx,0,sizeof(double)*nx*nmaxroot);
        }
    }
    
    /*
     * ----------------------------------------------------------------
     * PREPARE OUTPUT
     * ----------------------------------------------------------------
     */
    
    if(mxGetField(prhs[0], 0 ,"t")) { tdata = mxGetPr(mxGetField(prhs[0], 0 ,"t")); } else { mexErrMsgTxt("t not specified as field in solution struct!"); }
    
    if(mxGetField(prhs[0], 0 ,"x")) { xdata = mxGetPr(mxGetField(prhs[0], 0 ,"x")); } else { mexErrMsgTxt("x not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"y")) { ydata = mxGetPr(mxGetField(prhs[0], 0 ,"y")); } else { mexErrMsgTxt("y not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"xdot")) { xdotdata = mxGetPr(mxGetField(prhs[0], 0 ,"xdot")); } else { mexErrMsgTxt("xdot not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"J")) { Jdata = mxGetPr(mxGetField(prhs[0], 0 ,"J")); } else { mexErrMsgTxt("J not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"llh")) { llhdata = mxGetPr(mxGetField(prhs[0], 0 ,"llh")); } else { mexErrMsgTxt("llh not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"chi2")) { chi2data = mxGetPr(mxGetField(prhs[0], 0 ,"chi2")); } else { mexErrMsgTxt("chi2 not specified as field in solution struct!"); }
    if (nr>0) {
        if(mxGetField(prhs[0], 0 ,"root")) { rootdata = mxGetPr(mxGetField(prhs[0], 0 ,"root")); } else { mexErrMsgTxt("root not specified as field in solution struct!"); }
        if(mxGetField(prhs[0], 0 ,"rootval")) { rootvaldata = mxGetPr(mxGetField(prhs[0], 0 ,"rootval")); } else { mexErrMsgTxt("root not specified as field in solution struct!"); }
    }
    
    if(mxGetField(prhs[0], 0 ,"numsteps")) { numstepsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numsteps")); } else { mexErrMsgTxt("numsteps not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevals")) { numrhsevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevals")); } else { mexErrMsgTxt("numrhsevals not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numlinsolvsetups")) { numlinsolvsetupsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numlinsolvsetups")); } else { mexErrMsgTxt("numlinsolvsetups not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numerrtestfails")) { numerrtestfailsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numerrtestfails")); } else { mexErrMsgTxt("numerrtestfails not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"order")) { orderdata = mxGetPr(mxGetField(prhs[0], 0 ,"order")); } else { mexErrMsgTxt("order not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numnonlinsolviters")) { numnonlinsolvitersdata = mxGetPr(mxGetField(prhs[0], 0 ,"numnonlinsolviters")); } else { mexErrMsgTxt("numnonlinsolviters not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numjacevals")) { numjacevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numjacevals")); } else { mexErrMsgTxt("numjacevals not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numliniters")) { numlinitersdata = mxGetPr(mxGetField(prhs[0], 0 ,"numliniters")); } else { mexErrMsgTxt("numliniters not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numconvfails")) { numconvfailsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numconvfails")); } else { mexErrMsgTxt("numconvfails not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numprecevals")) { numprecevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numprecevals")); } else { mexErrMsgTxt("numprecevals not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numprecsolves")) { numprecsolvesdata = mxGetPr(mxGetField(prhs[0], 0 ,"numprecsolves")); } else { mexErrMsgTxt("numprecsolves not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numjtimesevals")) { numjtimesevalsdata = mxGetPr(mxGetField(prhs[0], 0 ,"numjtimesevals")); } else { mexErrMsgTxt("numjtimesevals not specified as field in solution struct!"); }
    
    if(mxGetField(prhs[0], 0 ,"numstepsS")) { numstepsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numstepsS")); } else { mexErrMsgTxt("numstepsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numrhsevalsS")) { numrhsevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numrhsevalsS")); } else { mexErrMsgTxt("numrhsevalsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numlinsolvsetupsS")) { numlinsolvsetupsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numlinsolvsetupsS")); } else { mexErrMsgTxt("numlinsolvsetupsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numerrtestfailsS")) { numerrtestfailsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numerrtestfailsS")); } else { mexErrMsgTxt("numerrtestfailsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"orderS")) { orderSdata = mxGetPr(mxGetField(prhs[0], 0 ,"orderS")); } else { mexErrMsgTxt("orderS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numnonlinsolvitersS")) { numnonlinsolvitersSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numnonlinsolvitersS")); } else { mexErrMsgTxt("numnonlinsolvitersS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numjacevalsS")) { numjacevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numjacevalsS")); } else { mexErrMsgTxt("numjacevalsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numlinitersS")) { numlinitersSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numlinitersS")); } else { mexErrMsgTxt("numlinitersS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numconvfailsS")) { numconvfailsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numconvfailsS")); } else { mexErrMsgTxt("numconvfailsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numprecevalsS")) { numprecevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numprecevalsS")); } else { mexErrMsgTxt("numprecevalsS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numprecsolvesS")) { numprecsolvesSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numprecsolvesS")); } else { mexErrMsgTxt("numprecsolvesS not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"numjtimesevalsS")) { numjtimesevalsSdata = mxGetPr(mxGetField(prhs[0], 0 ,"numjtimesevalsS")); } else { mexErrMsgTxt("numjtimesevalsS not specified as field in solution struct!"); }
    
    if (sensi >= 1) {
        if(mxGetField(prhs[0], 0 ,"llhS")) { llhSdata = mxGetPr(mxGetField(prhs[0], 0 ,"llhS")); } else { mexErrMsgTxt("llhS not specified as field in solution struct!"); }
        if (sensi >= 2) {
            if(mxGetField(prhs[0], 0 ,"llhS2")) { llhS2data = mxGetPr(mxGetField(prhs[0], 0 ,"llhS2")); } else { mexErrMsgTxt("llhS not specified as field in solution struct!"); }
        }
    }
    
    
    /*
     * ----------------------------------------------------------------
     * CVODE
     * ----------------------------------------------------------------
     */
    
    cv_status = 0.0;
    
    
    
    
    /* FORWARD PROBLEM */
    
    ncheck = 0; /* the number of (internal) checkpoints stored so far */
    
    t = tstart;
    
    nroots = 0;
    idisc = 0;
    for (it=0; it < nt; it++) {
        /* only integrate if no errors occured */
        if(cv_status == 0.0) {
            if(ts[it] > tstart) {
                if(nx>0) {
                    while (t<ts[it]) {
                        if (nr+ndisc>0) {
                            /* we have to find roots */
                            if(sensi >= 1) {
                                cv_status = CVodeF(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL, &ncheck);
                                if (cv_status==CV_ROOT_RETURN) {
                                    
                                    *status = CVodeGetRootInfo(cvode_mem, rootsfound);
                                    if (*status != CV_SUCCESS) goto error_return;
                                    /* EVENT ROOTS */
                                    if (nroots<nmaxroot) {
                                        for (ir=0; ir<nr; ir++){
                                            if(rootsfound[ir] == 1) {
                                                rootdata[nroots + nmaxroot*ir] = t;
                                                rootvaldata[nroots + nmaxroot*ir] = 0;
                                                rootidx[nroots] = ir;
                                                x_tmp = NV_DATA_S(x);
                                                /* extract sensitivity information */
                                                if(!mxIsNaN(mt[ir*nmaxroot+nroots])) {
                                                    if (event_model == LW_NORMAL) {
                                                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                                                        *chi2data += pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                                                    }
                                                    if (event_model == LW_NORMAL) {
                                                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+nroots] , 2);
                                                    }
                                                    if (sensi>=1) {
                                                        x_tmp = NV_DATA_S(x);
                                                        fdtdp(t,dtdp,data->p,data->k,x_tmp,plist,np,nmaxroot);
                                                        fdtdx(t,dtdx,data->p,data->k,x_tmp);
                                                        fdrvaldp(t,drvaldp,data->p,data->k,x_tmp,plist,np,nmaxroot);
                                                        fdrvaldx(t,drvaldx,data->p,data->k,x_tmp);
                                                        for (ip=0; ip<np; ip++) {
                                                            if(event_model == LW_NORMAL) {
                                                                drdp[ip] += ( dtdp[ip*nr+ir]* ( t - mt[ir*nmaxroot+nroots] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                                            }
                                                        }
                                                        for (ix=0; ix<nx; ix++) {
                                                            if(event_model  == LW_NORMAL) {
                                                                drdx[nroots+ix*nmaxroot] += ( dtdx[ix*nr+ir] * ( t - mt[ir*nmaxroot+nroots] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                                            }
                                                        }
                                                        for (ip=0; ip<np; ip++) {
                                                            if(event_model  == LW_NORMAL) {
                                                                drdp[ip] += ( drvaldp[ip*nr+ir]* ( 0 ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                                            }
                                                        }
                                                        for (ix=0; ix<nx; ix++) {
                                                            if(event_model  == LW_NORMAL) {
                                                                drdx[nroots+ix*nmaxroot] += ( drvaldx[ix*nr+ir] * ( 0 ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                                                            }
                                                        }
                                                    }
                                                }
                                                nroots++;
                                            }
                                        }
                                    }
                                    /* ROOTS FOR DISCONTINUITIES */
                                    for (ir=nr; ir<(nr+ndisc); ir++) {
                                        if(rootsfound[ir] != 0) {
                                            deltadisc(t,ir-nr,x,data);
                                            *status = CVodeReInit(cvode_mem,t,x);
                                            if (*status != CV_SUCCESS) goto error_return;
                                            if(idisc<nmaxdisc) {
                                                discs[idisc] = t;
                                                irdiscs[idisc] = ir-nr;
                                                idisc++;
                                            }
                                        }
                                    }
                                    /* if the root coincides with one of the output timepoints, we want to continue as if we reached that timepoint after accounting for the root */
                                    if (t==ts[it]) {
                                        cv_status = 0.0;
                                    }
                                }
                                if (cv_status == -22) {
                                    /* clustering of roots => turn off rootfinding */
                                    CVodeRootInit(cvode_mem, 0, NULL);
                                    cv_status = 0;
                                }
                            } else {
                                cv_status = CVode(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL);
                                if (cv_status==CV_ROOT_RETURN) {
                                    
                                    *status = CVodeGetRootInfo(cvode_mem, rootsfound);
                                    if (*status != CV_SUCCESS) goto error_return;
                                    if (nroots<nmaxroot) {
                                        for (ir=0; ir<nr; ir++){
                                            if(rootsfound[ir] == 1) {
                                                rootdata[nroots + nmaxroot*ir] = t;
                                                rootvaldata[nroots + nmaxroot*ir] = 0;
                                                rootidx[nroots] = ir;
                                            }
                                            if(!mxIsNaN(mt[ir*nmaxroot+nroots])) {
                                                if (event_model == LW_NORMAL) {
                                                    r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                                                    *chi2data += pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                                                }
                                                
                                                if (event_model == LW_NORMAL) {
                                                    r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( 0 )/tsigma[ir*nmaxroot+nroots] , 2);
                                                }
                                            }
                                            nroots++;
                                        }
                                    }
                                    
                                    /* ROOTS FOR DISCONTINUITIES */
                                    for (ir=nr; ir<(nr+ndisc); ir++) {
                                        if(rootsfound[ir] != 0) {
                                            deltadisc(t,ir-nr,x,data);
                                            *status = CVodeReInit(cvode_mem,t,x);
                                            if (*status != CV_SUCCESS) goto error_return;
                                            if(idisc<nmaxdisc) {
                                                discs[idisc] = t;
                                                irdiscs[idisc] = ir-nr;
                                                idisc++;
                                            }
                                        }
                                    }
                                    /* if the root coincides with one of the output timepoints, we want to continue as if we reached that timepoint after accounting for the root */
                                    if (t==ts[it]) {
                                        cv_status = 0.0;
                                    }
                                }
                                if (cv_status == -22) {
                                    /* clustering of roots => turn off rootfinding */
                                    CVodeRootInit(cvode_mem, 0, NULL);
                                    cv_status = 0;
                                }
                                
                            }
                        } else {
                            /* solve for forward problem and store checkpoints */
                            if(sensi >= 1) {
                                cv_status = CVodeF(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL, &ncheck);
                            } else {
                                cv_status = CVode(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL);
                            }
                        }
                        if (cv_status < 0) {
                            *status = cv_status;
                            if (*status != CV_SUCCESS) goto error_return;
                        }
                    }
                }
            }
            
            tdata[it] = ts[it];
            
            fxdot(t,x,xdot,data);
            
            x_tmp = NV_DATA_S(x);
            xdot_tmp = NV_DATA_S(xdot);
            /*fJ(nx,ts[it],x,xdot,Jtmp,data,NULL,NULL,NULL);*/
            for(ix=0; ix < nx; ix++) {
                xdata[ix*nt+it] = x_tmp[ix];
                xdotdata[ix*nt+it] = xdot_tmp[ix];
                /*for(jx=0; jx < nx;jx++){
                    Jdata[it+nt*(ix+nx*jx)] = Jtmp->data[ix+nx*jx];
                }*/
            }
            if(ts[it] > tstart) {
                *status = CVodeGetNumSteps(cvode_mem, &numsteps);
                if (*status != CV_SUCCESS) goto error_return;
                numstepsdata[it] = (double)numsteps;
                
                *status = CVodeGetNumRhsEvals(cvode_mem, &numrhsevals);
                if (*status != CV_SUCCESS) goto error_return;
                numrhsevalsdata[it] = (double)numrhsevals;
                
                *status = CVodeGetNumLinSolvSetups(cvode_mem, &numlinsolvsetups);
                if (*status != CV_SUCCESS) goto error_return;
                numlinsolvsetupsdata[it] = (double)numlinsolvsetups;
                
                *status = CVodeGetNumErrTestFails(cvode_mem, &numerrtestfails);
                if (*status != CV_SUCCESS) goto error_return;
                numerrtestfailsdata[it] = (double)numerrtestfails;
                
                *status = CVodeGetLastOrder(cvode_mem, &order);
                if (*status != CV_SUCCESS) goto error_return;
                orderdata[it] = (double)order;
                
                *status = CVodeGetNonlinSolvStats(cvode_mem, &numnonlinsolviters, &numconvfails);
                if (*status != CV_SUCCESS) goto error_return;
                numnonlinsolvitersdata[it] = (double)numnonlinsolviters;
                numconvfailsdata[it] = (double)numconvfails;
                
                if (linsol == 1) {
                    CVDlsGetNumJacEvals(cvode_mem, &numjacevals);
                }
                
                numjacevalsdata[it] = (double)numjacevals;
                
                if (linsol == 6 || linsol == 7 || linsol == 8) {
                    CVSpilsGetNumLinIters(cvode_mem, &numliniters);
                    numlinitersdata[it] = (double)numliniters;
                    
                    CVSpilsGetNumConvFails(cvode_mem, &numconvfails);
                    numconvfailsdata[it] = (double)numconvfails;
                    
                    CVSpilsGetNumPrecEvals(cvode_mem, &numprecevals);
                    numprecevalsdata[it] = (double)numprecevals;
                    
                    CVSpilsGetNumPrecSolves(cvode_mem, &numprecsolves);
                    numprecsolvesdata[it] = (double)numprecsolves;
                    
                    CVSpilsGetNumJtimesEvals(cvode_mem, &numjtimesevals);
                    numjtimesevalsdata[it] = (double)numjtimesevals;
                }
            }
            
        } else {
            for(ix=0; ix < nx; ix++) xdata[ix*nt+it] = 0.0;
        }
        
        if(cv_status == 0.0) {
            fy(ts[it],nt,it,ydata,data->p,data->k,data->u,xdata);
            for (iy=0; iy<ny; iy++) {
                if (data_model == LW_NORMAL) {
                    if(!mxIsNaN(my[iy*nt+it])){
                        g += 0.5*log(2*pi*pow(ysigma[iy*nt+it],2)) + 0.5*pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/ysigma[iy*nt+it] , 2);
                        *chi2data += pow( ( ydata[iy*nt+it] - my[iy*nt+it] )/ysigma[iy*nt+it] , 2);
                    }
                }
                if (data_model == LW_LOGNORMAL) {
                    if(!mxIsNaN(my[iy*nt+it])){
                        g += 0.5*log(2*pi*pow(ysigma[iy*nt+it]*ydata[iy*nt+it],2)) + 0.5*pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/ysigma[iy*nt+it] , 2);
                        *chi2data += pow( ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) )/ysigma[iy*nt+it] , 2);
                    }
                }
                if (data_model == LW_ONEOUTPUT) {
                    g += ydata[iy*nt+it];
                }
            }
            if (sensi >= 1) {
                for (ix=0; ix<nx; ix++) {
                    xnow[ix] = xdata[ix*nt + it];
                }
                for (iy=0; iy<ny; iy++) {
                    ynow[iy] = ydata[iy*nt + it];
                }
                fdydx(ts[it],dydx,ynow,data->p,data->k,xnow);
                fdydp(ts[it],dydp,ynow,data->p,data->k,xnow,plist,np,ny);
                for (iy=0; iy<ny; iy++) {
                    for (ip=0; ip<np; ip++) {
                        if(data_model == LW_NORMAL) {
                            if(!mxIsNaN(my[iy*nt+it])){
                                dgdp[ip] += ( dydp[ip*ny+iy]* ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( ysigma[iy*nt+it] , 2);
                            }
                        }
                        if(data_model == LW_LOGNORMAL) {
                            if(!mxIsNaN(my[iy*nt+it])){
                                dgdp[ip] += 1/(2*pi)*dydp[ip*ny+iy]/ydata[iy*nt+it] + ( dydp[ip*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( ysigma[iy*nt+it] , 2);
                            }
                        }
                        if(data_model == LW_ONEOUTPUT) {
                            dgdp[ip] += dydp[ip*ny+iy];
                        }
                    }
                    for (ix=0; ix<nx; ix++) {
                        if(data_model == LW_NORMAL) {
                            if(!mxIsNaN(my[iy*nt+it])){
                                dgdx[it+ix*nt] += ( dydx[ix*ny+iy] * ( ydata[iy*nt+it] - my[iy*nt+it] ) )/pow( ysigma[iy*nt+it] , 2);
                            }
                        }
                        if(data_model == LW_LOGNORMAL) {
                            if(!mxIsNaN(my[iy*nt+it])){
                                dgdx[it+ix*nt] += 1/(2*pi)*dydx[ix*ny+iy]/ydata[iy*nt+it] + ( dydx[ix*ny+iy]/ydata[iy*nt+it] * ( log(ydata[iy*nt+it]) - log(my[iy*nt+it]) ) )/pow( ysigma[iy*nt+it] , 2);
                            }
                        }
                        if(data_model == LW_ONEOUTPUT) {
                            dgdx[it+ix*nt] += dydx[ix*ny+iy];
                        }
                    }
                }
            }
        }
    }
    
    /* add event if we did not have one yet */
    if (nr>0) {
        if (nroots==0) {
            for (ir=0; ir<nr; ir++){
                rootdata[nroots + nmaxroot*ir] = t;
                froot(t, x, rootvaltmp, data);
                rootvaldata[nroots + nmaxroot*ir] = rootvaltmp[ir];
                rootidx[nroots] = ir;
                if(!mxIsNaN(mt[ir*nmaxroot+nroots])) {
                    if (event_model == LW_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                        *chi2data += pow( ( t - mt[ir*nmaxroot+nroots] )/tsigma[ir*nmaxroot+nroots] , 2);
                    }
                    
                    if (event_model == LW_NORMAL) {
                        r += 0.5*log(2*pi*pow(tsigma[ir*nmaxroot+nroots],2)) + 0.5*pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+nroots] , 2);
                        *chi2data += pow( ( rootvaltmp[ir] )/tsigma[ir*nmaxroot+nroots] , 2);
                    }
                    if (sensi>=1) {
                        x_tmp = NV_DATA_S(x);
                        fdtdp(t,dtdp,data->p,data->k,x_tmp,plist,np,nmaxroot);
                        fdtdx(t,dtdx,data->p,data->k,x_tmp);
                        fdrvaldp(t,drvaldp,data->p,data->k,x_tmp,plist,np,nmaxroot);
                        fdrvaldx(t,drvaldx,data->p,data->k,x_tmp);
                        for (ip=0; ip<np; ip++) {
                            if(event_model == LW_NORMAL) {
                                drdp[ip] += ( dtdp[ip*nr+ir]* ( t - mt[ir*nmaxroot+nroots] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                            if(event_model  == LW_LOGNORMAL) {
                                drdp[ip] += 1/(2*pi)*dtdp[ip*nr+ir]/t + ( dtdp[ip*nr+ir]/t * ( log(t) - log(mt[ir*nmaxroot+nroots]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(event_model  == LW_NORMAL) {
                                drdx[nroots+ix*nmaxroot] += ( dtdx[ix*nr+ir] * ( t - mt[ir*nmaxroot+nroots] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                            if(event_model  == LW_LOGNORMAL) {
                                drdx[nroots+ix*nmaxroot] += 1/(2*pi)*dtdx[ix*nr+ir]/t + ( dtdx[ix*nr+ir]/t * ( log(t) - log(mt[ir*nmaxroot+nroots]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                        }
                        for (ip=0; ip<np; ip++) {
                            if(event_model  == LW_NORMAL) {
                                drdp[ip] += ( drvaldp[ip*nr+ir]* ( rootvaltmp[ir] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                            if(event_model  == LW_LOGNORMAL) {
                                drdp[ip] += 1/(2*pi)*drvaldp[ip*nr+ir]/rootvaltmp[ir] + ( drvaldp[ip*nr+ir]/rootvaltmp[ir] * ( log(rootvaltmp[ir]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                        }
                        for (ix=0; ix<nx; ix++) {
                            if(event_model  == LW_NORMAL) {
                                drdx[nroots+ix*nmaxroot] += ( drvaldx[ix*nr+ir] * ( rootvaltmp[ir] ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                            if(event_model  == LW_LOGNORMAL) {
                                drdx[nroots+ix*nmaxroot] += 1/(2*pi)*drvaldx[ix*nr+ir]/rootvaltmp[ir] + ( drvaldx[ix*nr+ir]/rootvaltmp[ir] * ( log(rootvaltmp[ir]) ) )/pow( tsigma[ir*nmaxroot+nroots] , 2);
                            }
                        }
                    }
                }
                nroots++;
            }
        }
    }
    
    /* evaluate likelihood */
    
    *llhdata = - g - r;
    
    if (sensi >= 1) {
        /* only set output sensitivities if no errors occured */
        if(cv_status == 0.0) {
            
            xB = N_VNew_Serial(nx);
            
            xQB = N_VNew_Serial(np);
            
            /* BACKWARD PROBLEM */
            
            /* write initial conditions */
            if (xB == NULL) goto error_return;
            xB_tmp = NV_DATA_S(xB);
            memset(xB_tmp,0,sizeof(realtype)*nx);
            for (ix=0; ix<nx; ix++) {
                xB_tmp[ix] += dgdx[nt-1+ix*nt];
            }
            
            if (xQB == NULL) goto error_return;
            xQB_tmp = NV_DATA_S(xQB);
            memset(xQB_tmp,0,sizeof(realtype)*np);
            
            /* create backward problem */
            if (lmmB>2||lmmB<1) {
                mexErrMsgTxt("Illegal value for lmm!");
            }
            if (iterB>2||iterB<1) {
                mexErrMsgTxt("Illegal value for iter!");
            }
            *status = CVodeCreateB(cvode_mem, lmmB, iterB, &which);
            if (*status != CV_SUCCESS) goto error_return;
            
            
            /* allocate memory for the backward problem */
            *status = cvodewrap_binit(cvode_mem, which, xB, t);
            if (*status != CV_SUCCESS) goto error_return;
            
            /* specify integration tolerances for backward problem */
            *status = CVodeSStolerancesB(cvode_mem, which, RCONST(cvodes_rtol), RCONST(cvodes_atol));
            if(*status != CV_SUCCESS) goto error_return;
            
            /* Attach user data */
            *status = CVodeSetUserDataB(cvode_mem, which, data);
            if(*status != CV_SUCCESS) goto error_return;
            
            /* Number of maximal internal steps */
            *status = CVodeSetMaxNumStepsB(cvode_mem, which, cvodes_maxsteps);
            if(*status != CV_SUCCESS) goto error_return;
            
            switch (linsol) {
                    
                    /* DIRECT SOLVERS */
                    
                case LW_DENSE:
                    *status = CVDenseB(cvode_mem, which, nx);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                case LW_BAND:
                    *status = CVBandB(cvode_mem, which, nx, ubw, lbw);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetBandJacFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                case LW_LAPACKDENSE:
                    
                    /*#if SUNDIALS_BLAS_LAPACK
                     *status = CVLapackDenseB(cvode_mem, which, nx);
                     if (*status != CV_SUCCESS) goto error_return;
                     
                     *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                     if (*status != CV_SUCCESS) goto error_return;
                     #else*/
                    mexErrMsgTxt("Solver currently not supported!");
                    /*#endif*/
                    break;
                    
                case LW_LAPACKBAND:
                    
                    
                    /*#if SUNDIALS_BLAS_LAPACK
                     *status = CVLapackBandB(cvode_mem, which, nx, ubw, lbw);
                     if (*status != CV_SUCCESS) goto error_return;
                     
                     *status = cvodewrap_SetBandJacFnB(cvode_mem, which);
                     if (*status != CV_SUCCESS) goto error_return;
                     #else*/
                    mexErrMsgTxt("Solver currently not supported!");
                    /*#endif*/
                    break;
                    break;
                    
                case LW_DIAG:
                    *status = CVDiagB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                    /* ITERATIVE SOLVERS */
                    
                case LW_SPGMR:
                    *status = CVSpgmrB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                case LW_SPBCG:
                    *status = CVSpbcgB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                case LW_SPTFQMR:
                    *status = CVSptfqmrB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    break;
                    
                case LW_KLU:
                    *status = CVKLUB(cvode_mem, which, nx, nnz);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    *status = cvodewrap_SetSparseJacFnB(cvode_mem, which);
                    if (*status != CV_SUCCESS) goto error_return;
                    
                    /**status = CVKLUSetOrderingB(cvode_mem, which, 2);
                     if (*status != CV_SUCCESS) goto error_return;*/
                    
                    break;
                    
                default:
                    break;
            }
            
            /* Initialise quadrature calculation */
            *status = cvodewrap_qbinit(cvode_mem, which, xQB);
            if (*status != CV_SUCCESS) goto error_return;
            
            /* Enable Quadrature Error Control */
            *status = CVodeSetQuadErrConB(cvode_mem, which, TRUE);
            if (*status != CV_SUCCESS) goto error_return;
            
            *status = CVodeQuadSStolerancesB(cvode_mem, which, RCONST(cvodes_rtol), RCONST(cvodes_atol));
            if(*status != CV_SUCCESS) goto error_return;
            
            *status = CVodeSetStabLimDetB(cvode_mem, which, stldet); /* activates stability limit detection */
            if(*status != CV_SUCCESS) goto error_return;
            
            
            nroots--;
            idisc--;
            for (it=nt-1; it>0; it--) {
                /* enter while */
                rootflag = TRUE;
                while (rootflag || discflag) {
                    rootflag = FALSE;
                    discflag = FALSE;
                    if(idisc>-1) {
                        if(discs[idisc]>=ts[it-1]) {
                            if (nroots>-1) {
                                if (rootdata[nroots + nmaxroot*rootidx[nroots]] > discs[idisc]) {
                                    rootflag = TRUE;
                                } else {
                                    discflag = TRUE;
                                }
                            } else {
                                discflag = TRUE;
                            }
                        }
                    }
                    if (discflag == FALSE && rootflag == FALSE) {
                        if (nroots>-1) {
                            if (rootdata[nroots + nmaxroot*rootidx[nroots]] >= ts[it-1]) {
                                rootflag = TRUE;
                            }
                        }
                    }
                    
                    if (discflag) {
                        if (discs[idisc] == ts[it]) {
                            *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeGetB(cvode_mem, which, &t, xB);
                            if (*status != CV_SUCCESS) goto error_return;
                            ideltadisc(t,irdiscs[idisc],xB,xQB,data);
                            *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            Bdeltadisc(t,irdiscs[idisc],xB,data);
                            *status = CVodeReInitB(cvode_mem, which, t, xB);
                            if (*status != CV_SUCCESS) goto error_return;
                            t = discs[idisc];
                            idisc--;
                        } else {
                            cv_status = CVodeB(cvode_mem, discs[idisc], CV_NORMAL);
                            t = discs[idisc];
                            *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeGetB(cvode_mem, which, &t, xB);
                            if (*status != CV_SUCCESS) goto error_return;
                            ideltadisc(t,irdiscs[idisc],xB,xQB,data);
                            *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                            if (*status != CV_SUCCESS) goto error_return;
                            Bdeltadisc(t,irdiscs[idisc],xB,data);
                            *status = CVodeReInitB(cvode_mem, which, t, xB);
                            if (*status != CV_SUCCESS) goto error_return;
                            t = discs[idisc];
                            idisc--;
                        }
                    }
                    if (rootflag) {
                        *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        cv_status = CVodeB(cvode_mem, rootdata[nroots + nmaxroot*rootidx[nroots]], CV_NORMAL);
                        t = rootdata[nroots + nmaxroot*rootidx[nroots]];
                        *status = CVodeGetB(cvode_mem, which, &t, xB);
                        xB_tmp = NV_DATA_S(xB);
                        for (ix=0; ix<nx; ix++) {
                            xB_tmp[ix] += drdx[nroots + ix*nmaxroot];
                        }
                        *status = CVodeReInitB(cvode_mem, which, t, xB);
                        if (*status != CV_SUCCESS) goto error_return;

                        t = rootdata[nroots + nmaxroot*rootidx[nroots]];
                        nroots--;
                    }
                }
                if(cv_status == 0.0) {
                    if (nx>0) {
                        /* solve for backward problems */
                        if (ts[it-1] < t) {
                            cv_status = CVodeB(cvode_mem, RCONST(ts[it-1]), CV_NORMAL);
                        }
                        *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeGetB(cvode_mem, which, &t, xB);
                        if (*status != CV_SUCCESS) goto error_return;
                        
                        xB_tmp = NV_DATA_S(xB);
                        for (ix=0; ix<nx; ix++) {
                            xB_tmp[ix] += dgdx[it-1+ix*nt];
                        }
                        
                        t = ts[it-1];
                        *status = CVodeReInitB(cvode_mem, which, t, xB);
                        if (*status != CV_SUCCESS) goto error_return;
                    }
                }
                
            }
            
            for (ir=0; ir<nr; ir++) {
                if (nroots>-1) {
                    if (rootdata[nroots + nmaxroot*ir]< t && rootdata[nroots + nmaxroot*ir]>ts[it-1]) {
                        cv_status = CVodeB(cvode_mem, rootdata[nroots + nmaxroot*ir], CV_NORMAL);
                        nroots--;
                        *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeQuadReInitB(cvode_mem, which, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeGetB(cvode_mem, which, &t, xB);
                        xB_tmp = NV_DATA_S(xB);
                        for (ix=0; ix<nx; ix++) {
                            xB_tmp[ix] += drdx[it+ix*nt];
                        }
                        t = rootdata[nroots + nmaxroot*ir];
                        *status = CVodeReInitB(cvode_mem, which, t, xB);
                        if (*status != CV_SUCCESS) goto error_return;
                    }
                }
            }
            
            if (t>tstart) {
                if(cv_status == 0.0) {
                    if (nx>0) {
                        /* solve for backward problems */
                        cv_status = CVodeB(cvode_mem, tstart, CV_NORMAL);
                        
                        *status = CVodeGetQuadB(cvode_mem, which, &t, xQB);
                        if (*status != CV_SUCCESS) goto error_return;
                        *status = CVodeGetB(cvode_mem, which, &t, xB);
                        if (*status != CV_SUCCESS) goto error_return;
                    }
                }
            }
            
            /* evaluate initial values */
            sx = N_VCloneVectorArray_Serial(np,x);
            if (sx == NULL) goto error_return;
            fsx0(sx, data);
            
            if(cv_status == 0.0) {
                

                
                xB_tmp = NV_DATA_S(xB);
                
                for (ip=0; ip<np; ip++) {
                    llhS0[ip] = 0.0;
                    sx_tmp = NV_DATA_S(sx[ip]);
                    for (ix = 0; ix < nx; ix++) {
                        llhS0[ip] = llhS0[ip] + xB_tmp[ix] * sx_tmp[ix];
                    }
                }
                
                
                xQB_tmp = NV_DATA_S(xQB);
                
                for(ip=0; ip < np; ip++) {
                    llhSdata[ip] = - llhS0[ip] - dgdp[ip] - drdp[ip] - xQB_tmp[ip];
                }
                
                 cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, which);
                 
                 *status = CVodeGetNumSteps(cvode_memB, &numsteps);
                 if (*status != CV_SUCCESS) goto error_return;
                 numstepsSdata[it] = (double)numsteps;
                 
                 *status = CVodeGetNumRhsEvals(cvode_memB, &numrhsevals);
                 if (*status != CV_SUCCESS) goto error_return;
                 numrhsevalsSdata[it] = (double)numrhsevals;
                 
                 *status = CVodeGetNumLinSolvSetups(cvode_memB, &numlinsolvsetups);
                 if (*status != CV_SUCCESS) goto error_return;
                 numlinsolvsetupsSdata[it] = (double)numlinsolvsetups;
                 
                 *status = CVodeGetNumErrTestFails(cvode_memB, &numerrtestfails);
                 if (*status != CV_SUCCESS) goto error_return;
                 numerrtestfailsSdata[it] = (double)numerrtestfails;
                 
                 *status = CVodeGetLastOrder(cvode_memB, &order);
                 if (*status != CV_SUCCESS) goto error_return;
                 orderSdata[it] = (double)order;
                 
                 *status = CVodeGetNonlinSolvStats(cvode_memB, &numnonlinsolviters, &numconvfails);
                 if (*status != CV_SUCCESS) goto error_return;
                 numnonlinsolvitersSdata[it] = (double)numnonlinsolviters;
                 numconvfailsSdata[it] = (double)numconvfails;
                
                 if (linsol == 1) {
                    CVDlsGetNumJacEvals(cvode_memB, &numjacevals);
                    numjacevalsSdata[it] = (double)numjacevals;
                 }

                
                 if (linsol == 6 || linsol == 7 || linsol == 8) {
                    CVSpilsGetNumLinIters(cvode_memB, &numliniters);
                    numlinitersSdata[it] = (double)numliniters;
                     
                    CVSpilsGetNumConvFails(cvode_memB, &numconvfails);
                    numconvfailsSdata[it] = (double)numconvfails;
                     
                    CVSpilsGetNumPrecEvals(cvode_memB, &numprecevals);
                    numprecevalsSdata[it] = (double)numprecevals;
                     
                    CVSpilsGetNumPrecSolves(cvode_memB, &numprecsolves);
                    numprecsolvesSdata[it] = (double)numprecsolves;
                     
                    CVSpilsGetNumJtimesEvals(cvode_memB, &numjtimesevals);
                    numjtimesevalsSdata[it] = (double)numjtimesevals;
                 }
                
            } else {
                for(ip=0; ip < np; ip++) {
                    llhSdata[ip] = 0.0;
                }
            }
        } else {
            for(ip=0; ip < np; ip++) {
                llhSdata[ip] = 0.0;
            }
        }
    }
    
    
    
    /*
     * ----------------------------------------------------------------
     * OUTPUT
     * ----------------------------------------------------------------
     */
    
    *status = cv_status;
    
    /* Free memory */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(xdot);
        CVodeFree(&cvode_mem);
        DestroyMat(Jtmp);
        if (nr+ndisc>0) {
            if(rootsfound) mxFree(rootsfound);
        }
        if (nr>0) {
            if(rootvaltmp) mxFree(rootvaltmp);
            if(rootidx) mxFree(rootidx);
        }
        if (ndisc>0) {
            if(discs) mxFree(discs);
            if(irdiscs) mxFree(irdiscs);
        }
        
        if (sensi >= 1) {
            N_VDestroyVectorArray_Serial(sx,np);
            if(llhS0)     mxFree(llhS0);
            if(ynow)    mxFree(ynow);
            if(xnow)    mxFree(xnow);
            if(dydx)    mxFree(dydx);
            if(dydp)    mxFree(dydp);
            if(dgdp)    mxFree(dgdp);
            if(dgdx)    mxFree(dgdx);
            if(xB)      N_VDestroy_Serial(xB);
            if(xQB)     N_VDestroy_Serial(xQB);
        }
    }
    
    if(data)    mxFree(data);
    if(plist)   mxFree(plist);
    
    
    return;
    
error_return:
    return;
}