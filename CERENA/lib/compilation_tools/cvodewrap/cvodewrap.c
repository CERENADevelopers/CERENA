/*
 *  MATLAB usage: [status,t,x,y,sx,sy] = simulate_"model"(tout,theta,kappa,options)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <cvodes/cvodes.h>
/*#include <cvodes/cvodes_lapack.h>*/
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
/*#include <sundials/sundials_lapack.h> */ /* LAPACK linear solvers */
#include <klu.h>
#include <amd.h>
#include <colamd.h>
#include <btf.h>

/* user functions */
#include "cvodewrapfunctions.c"

/* Accessor macros */
#define Ith(v, i)     NV_Ith_S(v, i-1)        /* i-th vector component i=1..nx */

#define CW_DENSE       1
#define CW_BAND        2
#define CW_LAPACKDENSE 3
#define CW_LAPACKBAND  4
#define CW_DIAG        5
#define CW_SPGMR       6
#define CW_SPBCG       7
#define CW_SPTFQMR     8
#define CW_KLU         9

#define CW_FSA 1
#define CW_ASA 2




int sensi;
double cvodes_atol,cvodes_rtol;
int cvodes_maxsteps;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int cv_status;
    double *status;
    int np, nx, ny, nt, nr, nnz, ndisc; /* integers with problem dimensionality */
    int ip, ix, jx, iy, it, ir; /* integers for indexing in loops */
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
    N_Vector *xBt; /* state vector array for backwards problem starting at individual timepoints */
    N_Vector *xQBt; /* state vector array for backwards integration problem starting at individual timepoints */
    DlsMat Jtmp;
    
    double *ts, *tdata; /* timepoints */
    double *xdotdata; /* returned vector containing differential states*/
    double *dxdotdpdata;
    double *dydxdata;
    double *dydpdata;
    double *Jdata;
    double *rootdata, *rootSdata, *rootS2data;
    double *rootvaldata, *rootvalSdata, *rootvalS2data;
    double *rootvaltmp; /* returned vector containing states and sensitivities */
    double *xdata, *xSdata; /* returned vector containing states and sensitivities */
    double *ydata, *ySdata; /* returned vector containing outputs and sensitivities */
    
    int nmaxroot; /* maximum number of roots to record */
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
    
    
    double *yS0;
    realtype *x_tmp, *sx_tmp; /* temporary storage for states and sensitivities */
    realtype *xdot_tmp; /* temporary storage for states and sensitivities */
    realtype *xBt_tmp, *xQBt_tmp; /* temporary storage for backward states */
    double *dydx, *dydp, *xend, *yend; /* temporary storage for initiation of backward states */
    double *dxdotdp; /* temporary storage for initiation of backward states */
    double *p;  /* parameter vector */
    int    *plist;  /* index set of parameters for which to compute sensitivities */
    double *plistdata; /* input for plist */
    double *pbar;
    double *xbar;
    
    long int    Nd;     /* number of integration steps between two consecutive checkpoints. */
    
    int ism; /* a flag used to select the sensitivity solution method. Its value can be CV SIMULTANEOUS or CV STAGGERED */
    int sensi_meth; /* CW_FSA for forward sensitivity analysis, CW_ASA for adjoint sensitivity analysis */
    int linsol; /* linear solver specification */
    int interpType; /* CV_POLYNOMIAL or CV_HERMITE */
    int lmm; /* FORWARD PROBLEM: specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.  */
    int iter; /* FORWARD PROBLEM: specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.  */
    int lmmB; /* BACKWARD PROBLEM: specifies the linear multistep method and may be one of two possible values: CV ADAMS or CV BDF.  */
    int iterB; /* BACKWARD PROBLEM: specifies the type of nonlinear solver iteration and may be either CV NEWTON or CV FUNCTIONAL.  */
    int which = 0; /* BACKWARD PROBLEM: contains the identifier assigned by cvodes for the newly created back- ward problem. Any call to CVode*B functions requires such an identifier. */
    int *whichit; /* index set of identifiers, needed for later extraction */
    bool stldet; /* flag controlling stability limit detection */
    double stldetdata;
    bool error_corr = TRUE;
    int pretype= PREC_NONE; /* specifies the type of preconditioning and must be one of: PREC NONE, PREC LEFT, PREC RIGHT, or PREC BOTH. */
    double tlastroot; /* storage for last found root */
    
    double ss; /* flag for steady state analysis */
    
    bool b_sx0; /* flag which determines whether analytic sensitivities initialisation or provided initialisation should be used */
    double *sx0;
    
    
    
    /*
     * ----------------------------------------------------------------
     * EXTRACT INPUTS
     * ----------------------------------------------------------------
     */
    
    /* User data structure */
    data = (UserData) mxMalloc(sizeof *data);
    if (data == NULL) goto error_return;
    
    
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
    
    /* options */
    if (!prhs[4]) {
        mexErrMsgTxt("No options provided!");
    }
    
    if(mxGetField(prhs[4], 0 ,"tstart")) { tstart = mxGetScalar(mxGetField(prhs[4], 0 ,"tstart")); } else { mexErrMsgTxt("Parameter tstart not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"cvodes_atol")) { cvodes_atol = mxGetScalar(mxGetField(prhs[4], 0 ,"cvodes_atol")); } else { mexErrMsgTxt("Parameter cvodes_atol not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"cvodes_rtol")) { cvodes_rtol = mxGetScalar(mxGetField(prhs[4], 0 ,"cvodes_rtol")); } else { mexErrMsgTxt("Parameter cvodes_rtol not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"cvodes_maxsteps")) { cvodes_maxsteps = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"cvodes_maxsteps")); } else { mexErrMsgTxt("Parameter cvodes_maxsteps not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"nx")) { nx = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nx")); } else { mexErrMsgTxt("Parameter nx not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ny")) { ny = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ny")); } else { mexErrMsgTxt("Parameter ny not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"np")) { np = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"np")); } else { mexErrMsgTxt("Parameter np not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nr")) { nr = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nr")); } else { mexErrMsgTxt("Parameter nr not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ndisc")) { ndisc = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ndisc")); } else { mexErrMsgTxt("Parameter ndisc not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"nnz")) { nnz = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nnz")); } else { mexErrMsgTxt("Parameter nnz not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"lmm")) { lmm = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"lmm")); } else {  mexErrMsgTxt("Parameter lmm not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"iter")) { iter = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"iter")); } else { mexErrMsgTxt("Parameter iter not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"lmmB")) { lmmB = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"lmmB")); } else { mexErrMsgTxt("Parameter lmmB not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"iterB")) { iterB = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"iterB")); } else { mexErrMsgTxt("Parameter iterB not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"Nd"))  { Nd = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"Nd")); } else { mexErrMsgTxt("Parameter Nd not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"interpType"))  { interpType = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"interpType")); } else { mexErrMsgTxt("Parameter interpType not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"linsol")) { linsol = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"linsol")); } else { mexErrMsgTxt("Parameter linsol not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"stldet")) { stldetdata = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"stldet")); } else { mexErrMsgTxt("Parameter stldetdata not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"sensi")) { sensi = (int) mxGetScalar(mxGetField(prhs[4], 0 ,"sensi")); } else { mexErrMsgTxt("Parameter sensi not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"ism")) { ism = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ism")); } else { mexErrMsgTxt("Parameter ism not specified as field in options struct!"); }
    if(mxGetField(prhs[4], 0 ,"sensi_meth")) { sensi_meth = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"sensi_meth")); } else { mexErrMsgTxt("Parameter sensi_meth not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"nmaxroot")) { nmaxroot = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"nmaxroot")); } else { mexErrMsgTxt("Parameter nmaxroot not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"ubw")) { ubw = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ubw")); } else { mexErrMsgTxt("Parameter ubw not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"lbw")) { lbw = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"lbw")); } else { mexErrMsgTxt("Parameter lbw not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"ss")) { ss = (int)mxGetScalar(mxGetField(prhs[4], 0 ,"ss")); } else { mexErrMsgTxt("Parameter ss not specified as field in options struct!"); }
    
    if(mxGetField(prhs[4], 0 ,"sx0")) { sx0 = mxGetPr(mxGetField(prhs[4], 0 ,"sx0")); b_sx0 = TRUE;} else { b_sx0 = FALSE;}
    if (b_sx0) {
        /* check dimensions */
        if(mxGetN(mxGetField(prhs[4], 0 ,"sx0")) != np) { mexErrMsgTxt("Number of rows in sx0 field does not agree with number of model parameters!"); }
        if(mxGetM(mxGetField(prhs[4], 0 ,"sx0")) != nx) { mexErrMsgTxt("Number of columns in sx0 field does not agree with number of model states!"); }
    }
    
    
    if ((int)stldetdata>0.5) {
        stldet = TRUE;
    } else {
        stldet = FALSE;
    }
    
    
    data->np = &RCONST(np);
    data->ny = &RCONST(ny);
    if(mxGetField(prhs[4], 0 ,"qPositiveX")) { data->qpositivex = mxGetData(mxGetField(prhs[4], 0, "qPositiveX")); } else { mexErrMsgTxt("Parameter qPositiveX not specified as field in options struct!"); }
    
    /* plist */
    if (!prhs[5]) {
        mexErrMsgTxt("No parameter list provided!");
    }
    
    if(prhs[5]) {
        plistdata = mxGetPr(prhs[5]);
    }
    
    plist = mxMalloc(np*sizeof(int));
    for (ip=0; ip<np; ip++) {
        plist[ip] = (int)plistdata[ip];
    }
    
    data->plist = plist;
    
    /* pbar */
    if (!prhs[6]) {
        mexErrMsgTxt("No parameter scales provided!");
    }
    
    pbar = mxGetPr(prhs[6]);
    
    /* xscale */
    if (!prhs[7]) {
        mexErrMsgTxt("No state scales provided!");
    }
    
    xbar = mxGetPr(prhs[7]);
    
    
    
    
    /*
     * ----------------------------------------------------------------
     * CVODEINIT
     * ----------------------------------------------------------------
     */
    
    
    
    if (nx > 0) {
        
        /* write initial conditions */
        x = N_VNew_Serial(nx);
        xdot = N_VNew_Serial(nx);
        Jtmp = NewDenseMat(nx,nx);
        
        if(nr>0) rootvaltmp = mxMalloc(nr*sizeof(double));
        if(nr+ndisc>0) rootsfound = mxMalloc((nr+ndisc)*sizeof(double));
        
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
                
            case CW_DENSE:
                *status = CVDense(cvode_mem, nx);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetDenseJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case CW_BAND:
                *status = CVBand(cvode_mem, nx, ubw, lbw);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetBandJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case CW_LAPACKDENSE:
                mexErrMsgTxt("Solver currently not supported!");
                /**status = CVLapackDense(cvode_mem, nx);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 *status = cvodewrap_SetDenseJacFn(cvode_mem);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 break;*/
                
            case CW_LAPACKBAND:
                
                mexErrMsgTxt("Solver currently not supported!");
                /**status = CVLapackBand(cvode_mem, nx);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 *status = cvodewrap_SetBandJacFn(cvode_mem);
                 if (*status != CV_SUCCESS) goto error_return;
                 
                 break;*/
                
            case CW_DIAG:
                *status = CVDiag(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
                /* ITERATIVE SOLVERS */
                
            case CW_SPGMR:
                *status = CVSpgmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case CW_SPBCG:
                *status = CVSpbcg(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case CW_SPTFQMR:
                *status = CVSptfqmr(cvode_mem, PREC_NONE, CVSPILS_MAXL);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetJacTimesVecFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVSpilsSetPrecType(cvode_mem, pretype);
                if (*status != CV_SUCCESS) goto error_return;
                
                break;
                
            case CW_KLU:
                *status = CVKLU(cvode_mem, nx, nnz);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = cvodewrap_SetSparseJacFn(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVKLUSetOrdering(cvode_mem, 2); /* natural ordering */
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
        
        if (sensi_meth == CW_FSA) {
            
            if(nx>0) {
                
                /* Set sensitivity initial conditions */
                
                sx = N_VCloneVectorArray_Serial(np, x);
                if (sx == NULL) goto error_return;
                
                if(!b_sx0) {
                    fsx0(sx, data);
                } else {
                    for (ip=0; ip<np; ip++) {
                        sx_tmp = NV_DATA_S(sx[plist[ip]]);
                        for (ix=0; ix<nx; ix++) {
                            sx_tmp[ix] = sx0[ix + nx*plist[ip]];
                        }
                    }
                }
                
                
                /* Activate sensitivity calculations */
                
                *status = cvodewrap_SensInit1(cvode_mem, np, ism, sx);
                if (*status != CV_SUCCESS) goto error_return;
                
                /* Set sensitivity analysis optional inputs */
                *status = CVodeSetSensParams(cvode_mem, p, pbar, plist);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVodeSetSensErrCon(cvode_mem, error_corr);
                if (*status != CV_SUCCESS) goto error_return;
                
                *status = CVodeSensEEtolerances(cvode_mem);
                if (*status != CV_SUCCESS) goto error_return;
            }
        }
        
        if (sensi_meth == CW_ASA) {
            
            if(nx>0) {
                /* Allocate space for the adjoint computation */
                *status = CVodeAdjInit(cvode_mem, Nd, interpType);
                if (*status != CV_SUCCESS) goto error_return;
                
                whichit = mxMalloc(nt*sizeof(int));
                xend = mxMalloc(nx*sizeof(double));
                yend = mxMalloc(ny*sizeof(double));
                dydx = mxMalloc(ny*nx*sizeof(double));
                dydp = mxMalloc(ny*np*nt*sizeof(double));
                memset(dydp,0,sizeof(double)*ny*np*nt);
                yS0 = mxMalloc(np*ny*sizeof(double));
            }
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
    if(mxGetField(prhs[0], 0 ,"dxdotdp")) { dxdotdpdata = mxGetPr(mxGetField(prhs[0], 0 ,"dxdotdp")); } else { mexErrMsgTxt("dxdotdp not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"J")) { Jdata = mxGetPr(mxGetField(prhs[0], 0 ,"J")); } else { mexErrMsgTxt("J not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"dydx")) { dydxdata = mxGetPr(mxGetField(prhs[0], 0 ,"dydx")); } else { mexErrMsgTxt("dydx not specified as field in solution struct!"); }
    if(mxGetField(prhs[0], 0 ,"dydp")) { dydpdata = mxGetPr(mxGetField(prhs[0], 0 ,"dydp")); } else { mexErrMsgTxt("dydp not specified as field in solution struct!"); }
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
        if(mxGetField(prhs[0], 0 ,"yS")) { ySdata = mxGetPr(mxGetField(prhs[0], 0 ,"yS")); } else { mexErrMsgTxt("yS not specified as field in solution struct!"); }
        if (nr>0) {
            if(mxGetField(prhs[0], 0 ,"rootS")) { rootSdata = mxGetPr(mxGetField(prhs[0], 0 ,"rootS")); } else { mexErrMsgTxt("rootS not specified as field in solution struct!"); }
            if(mxGetField(prhs[0], 0 ,"rootvalS")) { rootvalSdata = mxGetPr(mxGetField(prhs[0], 0 ,"rootvalS")); } else { mexErrMsgTxt("rootvalS not specified as field in solution struct!"); }
            if (sensi >= 2) {
                if(mxGetField(prhs[0], 0 ,"rootS2")) { rootS2data = mxGetPr(mxGetField(prhs[0], 0 ,"rootS2")); } else { mexErrMsgTxt("rootS2 not specified as field in solution struct!"); }
                if(mxGetField(prhs[0], 0 ,"rootvalS2")) { rootvalS2data = mxGetPr(mxGetField(prhs[0], 0 ,"rootvalS2")); } else { mexErrMsgTxt("rootvalS2 not specified as field in solution struct!"); }
            }
        }
        if (sensi_meth == CW_FSA) {
            if(mxGetField(prhs[0], 0 ,"xS")) { xSdata = mxGetPr(mxGetField(prhs[0], 0 ,"xS")); } else { mexErrMsgTxt("xS not specified as field in solution struct!"); }
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
    tlastroot = 0;
    for (it=0; it < nt; it++) {
        /* only integrate if no errors occured */
        if(cv_status == 0.0) {
            if(ts[it] > tstart) {
                if(nx>0) {
                    while (t<ts[it]) {
                        if (nr+ndisc>0) {
                            /* we have to find roots */
                            if(sensi_meth == CW_ASA && sensi >= 1) {
                                cv_status = CVodeF(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL, &ncheck);
                                if (cv_status==CV_ROOT_RETURN) {
                                    rootdata[ir + nroots*nr] = t;
                                    ir++;
                                }
                            } else {
                                cv_status = CVode(cvode_mem, RCONST(ts[it]), x, &t, CV_NORMAL);
                                if (cv_status==CV_ROOT_RETURN) {
                                    if (t == tlastroot) {
                                        /* we are stuck in a root => turn off rootfinding */
                                        /* at some point we should find a more intelligent solution here, and turn on rootfinding again after some time */
                                        CVodeRootInit(cvode_mem, 0, NULL);
                                        cv_status = 0;
                                    }
                                    tlastroot = t;
                                    *status = CVodeGetRootInfo(cvode_mem, rootsfound);
                                    if (*status != CV_SUCCESS) goto error_return;
                                    /* ROOTS FOR ROOTFUNCTION */
                                    if (nroots<nmaxroot) {
                                        for (ir=0; ir<nr; ir++){ /* only look for roots of the rootfunction not discontinuities */
                                            if(rootsfound[ir] != 0) {
                                                rootdata[nroots + nmaxroot*ir] = t;
                                                rootvaldata[nroots + nmaxroot*ir] = 0;
                                                /* extract sensitivity information */
                                                if(sensi >= 1) {
                                                    if(sensi_meth == CW_FSA) {
                                                        *status = CVodeGetSens(cvode_mem, &t, sx);
                                                        if (*status != CV_SUCCESS) goto error_return;
                                                        for (ip=0; ip<np; ip++) {
                                                            rootSdata[nroots + nmaxroot*(ip*nr + ir)] = fsroot(t,plist[ip],ir,x,sx[plist[ip]],data);
                                                            rootvalSdata[nroots + nmaxroot*(ip*nr + ir)] = 0;
                                                            if (sensi >= 2) {
                                                                for (jp=0; jp<np; jp++) {
                                                                    rootS2data[nroots + nmaxroot*((np*ip+jp)*nr + ir)] = fs2root(t,plist[ip],plist[jp],ir,x,sx[plist[jp]],data);
                                                                    rootvalS2data[nroots + nmaxroot*((np*ip+jp)*nr + ir)] = 0;
                                                                }
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
                                            /* take care of deltas */
                                            /* sdeltadisc updates both x and sx, this needs to be done in combination for second order sensitivities */
                                            if(sensi >= 1) {
                                                if(sensi_meth == CW_FSA) {
                                                    *status = CVodeGetSens(cvode_mem, &t, sx);
                                                    if (*status != CV_SUCCESS) goto error_return;
                                                    sdeltadisc(t,ir-nr,x,sx,data);
                                                    CVodeSensReInit(cvode_mem, ism, sx);
                                                } else {
                                                    deltadisc(t,ir-nr,x,data);
                                                }
                                            } else {
                                                deltadisc(t,ir-nr,x,data);
                                            }
                                            /* reinitialize */
                                            *status = CVodeReInit(cvode_mem,t,x);
                                            if (*status != CV_SUCCESS) goto error_return;
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
                            if(sensi_meth == CW_ASA && sensi >= 1) {
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
            x_tmp = NV_DATA_S(x);
            for (ix=0; ix<nx; ix++) {
                xdata[it+nt*ix] = x_tmp[ix];
                if(data->qpositivex[ix]>0.5 && xdata[ix*nt+it]<0.0) xdata[ix*nt+it] = 0.0;
            }
            
            if (it == nt-1) {
                if( ss > 0) {
                    fxdot(t,x,xdot,data);
                    xdot_tmp = NV_DATA_S(xdot);
                    fJ(nx,ts[it],x,xdot,Jtmp,data,NULL,NULL,NULL);
                    for(ix=0; ix < nx; ix++) {
                        xdotdata[ix] = xdot_tmp[ix];
                        for(jx=0; jx < nx;jx++){
                            Jdata[ix+nx*jx] = Jtmp->data[ix+nx*jx];
                        }
                        /* set negative values to zeros */
                    }
                    fdxdotdp(t,nt,it,dxdotdpdata,data->p,data->k,data->u,xdata,plist,np,nx);
                    fdydx(ts[it],dydxdata,yend,data->p,data->k,xend);
                    fdydp(ts[it],nt,0,dydpdata,ydata,data->p,data->k,data->u,xdata,plist,np,ny);
                }
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
                	*status = CVDlsGetNumJacEvals(cvode_mem, &numjacevals);
                	if (*status != CV_SUCCESS) goto error_return;
                	numjacevalsdata[it] = (double)numjacevals;
                
                    *status = CVSpilsGetNumLinIters(cvode_mem, &numliniters);
                    if (*status != CV_SUCCESS) goto error_return;
                    numlinitersdata[it] = (double)numliniters;
                }
                
                if (linsol == 6 || linsol == 7 || linsol == 8) {
                    *status = CVSpilsGetNumConvFails(cvode_mem, &numconvfails);
                    if (*status != CV_SUCCESS) goto error_return;
                    numconvfailsdata[it] = (double)numconvfails;
                    
                    *status = CVSpilsGetNumPrecEvals(cvode_mem, &numprecevals);
                    if (*status != CV_SUCCESS) goto error_return;
                    numprecevalsdata[it] = (double)numprecevals;
                    
                    *status = CVSpilsGetNumPrecSolves(cvode_mem, &numprecsolves);
                    if (*status != CV_SUCCESS) goto error_return;
                    numprecsolvesdata[it] = (double)numprecsolves;
                    
                    *status = CVSpilsGetNumJtimesEvals(cvode_mem, &numjtimesevals);
                    if (*status != CV_SUCCESS) goto error_return;
                    numjtimesevalsdata[it] = (double)numjtimesevals;
                }
                
            }
            
        } else {
            for(ix=0; ix < nx; ix++) xdata[ix*nt+it] = 0.0;
        }
        
        if(cv_status == 0.0) {
            fy(ts[it],nt,it,ydata,data->p,data->k,data->u,xdata);
            if (sensi >= 1) {
                if (sensi_meth == CW_ASA) {
                    fdydp(ts[it],nt,it,dydp,ydata,data->p,data->k,data->u,xdata,plist,np,ny);
                }
                if (sensi_meth == CW_FSA) {
                    for(ip=0; ip < np; ip++) {
                        if(nx>0) {
                            if(ts[it] > tstart) {
                                *status = CVodeGetSens(cvode_mem, &t, sx);
                                if (*status != CV_SUCCESS) goto error_return;
                            }
                            
                            sx_tmp = NV_DATA_S(sx[ip]);
                            for(ix=0; ix < nx; ix++) {
                                xSdata[(ip*nx + ix)*nt + it] = sx_tmp[ix];
                            }
                        }
                    }
                    
                    *status = CVodeGetSensNumRhsEvals(cvode_mem, &numrhsevals);
                    if (*status != CV_SUCCESS) goto error_return;
                    numrhsevalsSdata[it] = (double)numrhsevals;
                    
                    *status = CVodeGetSensNumLinSolvSetups(cvode_mem, &numlinsolvsetups);
                    if (*status != CV_SUCCESS) goto error_return;
                    numlinsolvsetupsSdata[it] = (double)numlinsolvsetups;
                    
                    *status = CVodeGetSensNumErrTestFails(cvode_mem, &numerrtestfails);
                    if (*status != CV_SUCCESS) goto error_return;
                    numerrtestfailsSdata[it] = (double)numerrtestfails;
                    
                    
                    fsy(ts[it],nt,it,plist,nx,ny,np,ySdata,data->p,data->k,xdata,xSdata);
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
                /* extract sensitivity information */
                if(sensi >= 1) {
                    if(sensi_meth == CW_FSA) {
                        *status = CVodeGetSens(cvode_mem, &t, sx);
                        if (*status != CV_SUCCESS) goto error_return;
                        for (ip=0; ip<np; ip++) {
                            rootSdata[nroots + nmaxroot*(ip*nr + ir)] = 0;
                            rootvalSdata[nroots + nmaxroot*(ip*nr + ir)] = fsrootval(t,plist[ip],ir,x,sx[plist[ip]],data);
                            if (sensi >= 2) {
                                for (jp=0; jp<np; jp++) {
                                    rootS2data[nroots + nmaxroot*((np*ip+jp)*nr + ir)] = 0;
                                    rootvalS2data[nroots + nmaxroot*((np*ip+jp)*nr + ir)] = fs2rootval(t,plist[ip],plist[jp],ir,x,sx[plist[jp]],data);
                                }
                            }
                        }
                    }
                }
            }
            nroots++;
        }
    }
    
    
    if (sensi >= 1) {
        if (sensi_meth == CW_ASA) {
            /* only set output sensitivities if no errors occured */
            if(cv_status == 0.0) {
                
                xB = N_VNew_Serial(ny*nx);
                xBt = N_VCloneVectorArray_Serial(nt,xB);
                
                xQB = N_VNew_Serial(ny*np);
                xQBt = N_VCloneVectorArray_Serial(nt,xQB);
                
                for (it=0; it < nt; it++) {
                    
                    t = ts[it];
                    
                    /* BACKWARD PROBLEM */
                    
                    for (ix=0; ix<nx; ix++) {
                        xend[ix] = xdata[ix*nt + it];
                    }
                    for (iy=0; iy<ny; iy++) {
                        yend[iy] = ydata[iy*nt + it];
                    }
                    
                    fdydx(ts[it],dydx,yend,data->p,data->k,xend);
                    
                    /* write initial conditions */
                    if (xBt[it] == NULL) goto error_return;
                    xBt_tmp = NV_DATA_S(xBt[it]);
                    for (ix=0; ix<nx; ix++) {
                        for (iy=0; iy<ny; iy++) {
                            xBt_tmp[ix*ny+iy] = dydx[ix*ny+iy];
                        }
                    }
                    
                    if (xQBt[it] == NULL) goto error_return;
                    xQBt_tmp = NV_DATA_S(xQBt[it]);
                    memset(xQBt_tmp,0,sizeof(realtype)*ny*np);
                    
                    if(it>tstart){
                        /* create backward problem */
                        if (lmmB>2||lmmB<1) {
                            mexErrMsgTxt("Illegal value for lmm!");
                        }
                        if (iterB>2||iterB<1) {
                            mexErrMsgTxt("Illegal value for iter!");
                        }
                        *status = CVodeCreateB(cvode_mem, lmmB, iterB, &which);
                        if (*status != CV_SUCCESS) goto error_return;
                        
                        whichit[it] = which;
                        
                        /* allocate memory for the backward problem */
                        *status = cvodewrap_binit(cvode_mem, which, xBt[it], t);
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
                                
                            case CW_DENSE:
                                *status = CVDenseB(cvode_mem, which, nx*ny);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_BAND:
                                /**status = CVBandB(cvode_mem, which, nx*ny);*/
                                *status = -1;
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetBandJacFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_LAPACKDENSE:
                                mexErrMsgTxt("Solver currently not supported!");
                                /**status = CVLapackDenseB(cvode_mem, which, nx*ny);
                                 if (*status != CV_SUCCESS) goto error_return;
                                 
                                 *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                                 if (*status != CV_SUCCESS) goto error_return;
                                 
                                 break;*/
                                
                            case CW_LAPACKBAND:
                                /**status = CVLapackBandB(cvode_mem, which, nx*ny);*/
                                *status = -1;
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetBandJacFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_DIAG:
                                /**status = CVDiagB(cvode_mem, which);*/
                                *status = -1;
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetDenseJacFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                                /* ITERATIVE SOLVERS */
                                
                            case CW_SPGMR:
                                *status = CVSpgmrB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_SPBCG:
                                *status = CVSpbcgB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_SPTFQMR:
                                *status = CVSptfqmrB(cvode_mem, which, PREC_NONE, CVSPILS_MAXL);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetJacTimesVecFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = CVSpilsSetPrecTypeB(cvode_mem, which, pretype);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            case CW_KLU:
                                *status = CVKLUB(cvode_mem, which, nx, nnz);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                *status = cvodewrap_SetSparseJacFnB(cvode_mem, which);
                                if (*status != CV_SUCCESS) goto error_return;
                                
                                break;
                                
                            default:
                                break;
                        }
                        
                        
                        /* Initialise quadrature calculation */
                        *status = cvodewrap_qbinit(cvode_mem, which, xQBt[it]);
                        if (*status != CV_SUCCESS) goto error_return;
                        
                        *status = CVodeSetStabLimDetB(cvode_mem, which, stldet); /* activates stability limit detection */
                        if(*status != CV_SUCCESS) goto error_return;
                        
                        
                        
                    } else {
                        whichit[it] = -1;
                    }
                }
                if(cv_status == 0.0) {
                    if (ny*nx>0) {
                        /* solve for backward problems */
                        cv_status = CVodeB(cvode_mem, RCONST(tstart), CV_NORMAL);
                    }
                }
                
                /* evaluate initial values */
                sx = N_VCloneVectorArray_Serial(np,x);
                if (sx == NULL) goto error_return;
                fsx0(sx, data);
                
                for (it=0; it < nt; it++) {
                    if(cv_status == 0.0) {
                        
                        if (ts[it]>tstart) {
                            *status = CVodeGetQuadB(cvode_mem, whichit[it], &t, xQBt[it]);
                            if (*status != CV_SUCCESS) goto error_return;
                            *status = CVodeGetB(cvode_mem, whichit[it], &t, xBt[it]);
                            if (*status != CV_SUCCESS) goto error_return;
                        }
                        
                        xBt_tmp = NV_DATA_S(xBt[it]);
                        
                        for (ip=0; ip<np; ip++) {
                            for (iy = 0; iy<ny; iy++) {
                                yS0[iy + ip*ny] = 0.0;
                                sx_tmp = NV_DATA_S(sx[ip]);
                                for (ix = 0; ix < nx; ix++) {
                                    yS0[iy + ip*ny] = yS0[iy + ip*ny] + xBt_tmp[ix + iy*nx] * sx_tmp[ix];
                                }
                            }
                        }
                        
                        xQBt_tmp = NV_DATA_S(xQBt[it]);
                        
                        for(iy=0; iy < ny; iy++) {
                            for(ip=0; ip < np; ip++) {
                                ySdata[it + nt*(ip*ny + iy)] = yS0[iy + ip*ny] + dydp[it + nt*(iy + ip*ny)] + xQBt_tmp[iy + ip*ny];
                            }
                        }
                        
                        if (ts[it]>tstart) {
                            cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, whichit[it]);
                            
                            *status = CVodeGetNumSteps(cvode_memB, &numsteps);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { numstepsSdata[it] = numsteps; } else { numstepsSdata[it] = numstepsSdata[it-1] + (double)numsteps; }
                            
                            *status = CVodeGetNumRhsEvals(cvode_memB, &numrhsevals);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { numrhsevalsSdata[it] = (double)numrhsevals; } else { numrhsevalsSdata[it] = numrhsevalsSdata[it-1] + (double)numrhsevals; }
                            
                            *status = CVodeGetNumLinSolvSetups(cvode_memB, &numlinsolvsetups);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { numlinsolvsetupsSdata[it] = (double)numlinsolvsetups; } else { numlinsolvsetupsSdata[it] = numlinsolvsetupsSdata[it-1] + (double)numlinsolvsetups; }
                            
                            *status = CVodeGetNumErrTestFails(cvode_memB, &numerrtestfails);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { numerrtestfailsSdata[it] = (double)numerrtestfails; } else { numerrtestfailsSdata[it] = numerrtestfailsSdata[it-1] + (double)numerrtestfails; }
                            
                            *status = CVodeGetLastOrder(cvode_memB, &order);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { orderSdata[it] = (double)order; } else { orderSdata[it] = (double)order; }
                            
                            *status = CVodeGetNonlinSolvStats(cvode_memB, &numnonlinsolviters, &numconvfails);
                            if (*status != CV_SUCCESS) goto error_return;
                            if(it==0) { numnonlinsolvitersSdata[it] = (double)numnonlinsolviters; } else { numnonlinsolvitersSdata[it] = numnonlinsolvitersSdata[it-1] + (double)numnonlinsolviters; }
                            if(it==0) { numconvfailsSdata[it] = (double)numconvfails; } else { numconvfailsSdata[it] = numconvfailsSdata[it-1] + (double)numconvfails; }
                            
                            if (linsol == 1) {
                                *status = CVDlsGetNumJacEvals(cvode_memB, &numjacevals);
                                if(it==0) { numjacevalsSdata[it] = (double)numjacevals; } else { numjacevalsSdata[it] = numjacevalsSdata[it-1] + (double)numjacevals; }
                            }
                            
                            if (linsol == 6 || linsol == 7 || linsol == 8 ) {
                                *status = CVSpilsGetNumLinIters(cvode_memB, &numliniters);
                                if (*status != CV_SUCCESS) goto error_return;
                                if(it==0) { numlinitersSdata[it] = (double)numliniters; } else { numlinitersSdata[it] = numlinitersSdata[it-1] + (double)numliniters; }
                                
                                *status = CVSpilsGetNumConvFails(cvode_memB, &numconvfails);
                                if (*status != CV_SUCCESS) goto error_return;
                                if(it==0) { numconvfailsSdata[it] = (double)numconvfails; } else { numconvfailsSdata[it] = numconvfailsSdata[it-1] + (double)numconvfails; }
                                
                                *status = CVSpilsGetNumPrecEvals(cvode_memB, &numprecevals);
                                if (*status != CV_SUCCESS) goto error_return;
                                if(it==0) { numprecevalsSdata[it] = (double)numprecevals; } else { numprecevalsSdata[it] = numprecevalsSdata[it-1] + (double)numprecevals; }
                                
                                *status = CVSpilsGetNumPrecSolves(cvode_memB, &numprecsolves);
                                if (*status != CV_SUCCESS) goto error_return;
                                if(it==0) { numprecsolvesSdata[it] = (double)numprecsolves; } else { numprecsolvesSdata[it] = numprecsolvesSdata[it-1] + (double)numprecsolves; }
                                
                                *status = CVSpilsGetNumJtimesEvals(cvode_memB, &numjtimesevals);
                                if (*status != CV_SUCCESS) goto error_return;
                                if(it==0) { numjtimesevalsSdata[it] = (double)numjtimesevals; } else { numjtimesevalsSdata[it] = numjtimesevalsSdata[it-1] + (double)numjtimesevals; }
                            }
                            
                        }
                        
                    } else {
                        for(iy=0; iy < ny; iy++) {
                            for(ip=0; ip < np; ip++) {
                                ySdata[it + nt*(ip*ny + iy)] = 0.0;
                            }
                        }
                    }
                }
            } else {
                for (it=0; it < nt; it++) {
                    for(iy=0; iy < ny; iy++) {
                        for(ip=0; ip < np; ip++) {
                            ySdata[it + nt*(ip*ny + iy)] = 0.0;
                        }
                    }
                }
            }
        }
    }
    
    
    /*
     * ----------------------------------------------------------------
     * OUTPUT
     * ----------------------------------------------------------------
     */
    
    /*
     *  MATLAB usage: [status,t,x,y,sy] = simu_ode(tout,theta)
     */
    
    *status = cv_status;
    
    /* Free memory */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(xdot);
        CVodeFree(&cvode_mem);
        if (nr+ndisc>0) {
            if(rootsfound) mxFree(rootsfound);
        }
        if (nr>0) {
            if(rootvaltmp) mxFree(rootvaltmp);
        }
        if (sensi >= 1) {
            N_VDestroyVectorArray_Serial(sx,np);
            if (sensi_meth == CW_ASA) {
                if(yS0)     mxFree(yS0);
                if(yend)    mxFree(yend);
                if(xend)    mxFree(xend);
                if(dydx)    mxFree(dydx);
                if(dydp)    mxFree(dydp);
                if(whichit) mxFree(whichit);
                if(xB)      N_VDestroy_Serial(xB);
                if(xQB)     N_VDestroy_Serial(xQB);
                if(xBt)     N_VDestroyVectorArray_Serial(xBt,nt);
                if(xQBt)    N_VDestroyVectorArray_Serial(xQBt,nt);
            }
            
        }
    }
    
    if(data)    mxFree(data);
    if(plist)   mxFree(plist);
    
    
    return;
    
error_return:
    return;
}