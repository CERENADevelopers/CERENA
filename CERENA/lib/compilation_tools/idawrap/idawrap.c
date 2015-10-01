/*
 *  MATLAB usage: [status,t,x,y,sx,sy] = simu_ode(tout,theta,kappa,options)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <idas/idas.h>           /* prototypes for IDAS fcts. and consts. */
#include <idas/idas_dense.h>     /* prototype for IDASDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

/* Accessor macros */
#define Ith(v, i)     NV_Ith_S(v, i-1)        /* i-th vector component i=1..nx */
#define IJth(A, i, j) DENSE_ELEM(A, i-1, j-1) /* (i,j)-th matrix component i,j=1..nx */

/* user functions */
#include "idawrapfunctions.c"



int sensi;
double idas_atol,idas_rtol;
int idas_maxsteps;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int status, ida_status;
    int nx, ny, np, nt;
    int ix, iy, ip, it;
    
    int dims_xS[3], dims_yS[3];
    
    void *ida_mem;
    UserData data;
    
    realtype t;
    double tstart;
    N_Vector x, dx;
    N_Vector atols_ss;
    N_Vector *sx, *sdx;
    N_Vector id;
    realtype *sx_tmp, *x_tmp;
    realtype *sdx_tmp, *dx_tmp;
    realtype *id_tmp;
    
    double *qpositivex;
    
    double *ts, *tdata;
    double *xdata, *xSdata;
    double *dxdata, *dxSdata;
    double *ydata, *ySdata;
    double *dplist, *idlist;
    double *p;
    int    *plist;
    
    
    int sensi_meth = IDA_SIMULTANEOUS; /* IDA_SIMULTANEOUS or IDA_STAGGERED */
    bool error_corr = TRUE;
    
    
    idas_atol = mxGetScalar(mxGetField(prhs[3], 0 ,"idas_atol"));
    idas_rtol = mxGetScalar(mxGetField(prhs[3], 0 ,"idas_rtol"));
    idas_maxsteps = mxGetScalar(mxGetField(prhs[3], 0 ,"idas_maxsteps"));
    dplist = mxGetPr(prhs[4]);
    nx = mxGetNumberOfElements(mxGetField(prhs[3], 0, "xs"));
    ny = mxGetNumberOfElements(mxGetField(prhs[3], 0, "ys"));
    np = mxGetNumberOfElements(mxGetField(prhs[3], 0, "ps"));
    
    plist = mxMalloc(np*sizeof(int));
    
    for (ip=0; ip<np; ip++) {
        plist[ip] = (int)dplist[ip];
    }
    
    if (nlhs > 5) { /* sensitivities are on */
        sensi = 1;
    } else  { /* sensitivities are off */
        sensi = 0;
    }
    
    /*
     * ----------------------------------------------------------------
     * IDAINIT
     * ----------------------------------------------------------------
     */
    
    
    /* extract input time*/
    
    ts = mxGetPr(prhs[0]);
    tstart = mxGetScalar(mxGetField(prhs[3], 0 ,"tstart"));
    nt = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    
    /* User data structure */
    data = (UserData) mxMalloc(sizeof *data);
    if (data == NULL) goto error_return;
    
    data->qpositivex = mxGetData(mxGetField(prhs[3], 0, "qPositiveX"));
    data->u = mxGetData(mxGetField(prhs[3], 0, "uNum"));
    
    p = mxMalloc(np*sizeof(double));
    
    data->par = mxGetPr(prhs[1]);
    
    for (ip=0; ip<np; ip++) {
        p[ip] = data->par[plist[ip]];
    }
    
    data->p = p;
    data->plist = plist;
    data->k = mxGetPr(prhs[2]);
    
    data->v = mxGetData(mxGetField(prhs[3], 0, "vNum"));
    data->dvdx = mxGetData(mxGetField(prhs[3], 0, "dvdxNum"));
    data->dvddx = mxGetData(mxGetField(prhs[3], 0, "dvddxNum"));
    data->dvdu = mxGetData(mxGetField(prhs[3], 0, "dvduNum"));
    data->dvdp = mxGetData(mxGetField(prhs[3], 0, "dvdpNum"));
    
    /* fill for t=0 */
    fu(data, tstart);
    
    if (nx > 0) {
        
        /* write initial conditions */
        x = N_VNew_Serial(nx);
        x_tmp = NV_DATA_S(x);
        if (x == NULL) goto error_return;
        memset(x_tmp,0,sizeof(realtype)*nx);
        fx0(x, data);
        
        dx = N_VNew_Serial(nx);
        dx_tmp = NV_DATA_S(dx);
        if (dx == NULL) goto error_return;
        memset(dx_tmp,0,sizeof(realtype)*nx);
        fdx0(x, dx, data);
        
        fv(data, tstart, x, dx);
        
        /* Create IDAS object */
        ida_mem = IDACreate();
        if (ida_mem == NULL) goto error_return;
        
        /* Initialise IDAS solver */
        status = idawrap_init(ida_mem, x, dx, tstart);
        if (status != IDA_SUCCESS) goto error_return;
        
        /* Number of maximal internal steps */
        status = IDASetMaxNumSteps(ida_mem, idas_maxsteps);
        if(status != IDA_SUCCESS) goto error_return;
        
        /* Use private function to compute error weights */
        status = IDASStolerances(ida_mem, RCONST(idas_rtol), RCONST(idas_atol));
        if(status != IDA_SUCCESS) goto error_return;
        
        /* Attach user data */
        status = IDASetUserData(ida_mem, data);
        if(status != IDA_SUCCESS) goto error_return;
        
        /* Attach linear solver */
        status = IDADense(ida_mem, nx);
        if (status != IDA_SUCCESS) goto error_return;
        
        status = idawrap_SetDenseJacFn(ida_mem);
        if (status != IDA_SUCCESS) goto error_return;
    }
    
    /*
     * ----------------------------------------------------------------
     * IDASSENSINIT
     * ----------------------------------------------------------------
     */
    
    if ( sensi == 1) {
        
        data->su = mxGetData(mxGetField(prhs[3], 0, "suNum"));
        data->sv = mxGetData(mxGetField(prhs[3], 0, "svNum"));
        
        fsu(data, tstart);
        
        if(nx>0) {
            /* Load sensitivity initial conditions */
            sx = N_VCloneVectorArray_Serial(np, x);
            if (sx == NULL) goto error_return;
            for(ip=0; ip < np; ip++) {
                sx_tmp = NV_DATA_S(sx[ip]);
                memset(sx_tmp,0,sizeof(realtype)*nx);
            }
            for (ip=0;ip<np;ip++) fsx0(ip, sx[ip], x, dx, data);
            
            sdx = N_VCloneVectorArray_Serial(np, dx);
            if (sdx == NULL) goto error_return;
            for(ip=0; ip < np; ip++) {
                sdx_tmp = NV_DATA_S(sdx[ip]);
                memset(sdx_tmp,0,sizeof(realtype)*nx);
            }
            for (ip=0;ip<np;ip++) fsdx0(ip, sdx[ip], x, dx, data);
            
            
            fsv(data, tstart, x, dx);
            
            status = idawrap_SensInit(ida_mem, np, sensi_meth, sx, sdx);
            if (status != IDA_SUCCESS) goto error_return;
            
            atols_ss = N_VNew_Serial(np);
            if (atols_ss == NULL) {return;}
            for (ip=0; ip<np; ip++) Ith(atols_ss, ip+1) = idas_atol;
            
            status = IDASensSStolerances(ida_mem, RCONST(idas_rtol), N_VGetArrayPointer(atols_ss));
            if (status != IDA_SUCCESS) goto error_return;
            
            status = IDASetSensErrCon(ida_mem, error_corr);
            if (status != IDA_SUCCESS) goto error_return;
        }
        
    }
    
    /*
     * ----------------------------------------------------------------
     * PREPARE OUTPUT
     * ----------------------------------------------------------------
     */
    
    /*
     *  MATLAB usage: [status,t,x,dx,y,sx,sdx,sy] = simulate_modelname(tout,theta)
     */
    
    plhs[1] = mxCreateDoubleMatrix(1,nt,mxREAL);
    tdata = mxGetPr(plhs[1]);
    
    /* Solution vector(s) */
    
    plhs[2] = mxCreateDoubleMatrix(nt,nx,mxREAL);
    xdata = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(nt,nx,mxREAL);
    dxdata = mxGetPr(plhs[3]);
    
    plhs[4] = mxCreateDoubleMatrix(nt,ny,mxREAL);
    ydata = mxGetPr(plhs[4]);
    
    if (sensi == 1) {
        dims_xS[0] = nt;
        dims_xS[1] = nx;
        dims_xS[2] = np;
        plhs[5] = mxCreateNumericArray(3, dims_xS, mxDOUBLE_CLASS, mxREAL);
        xSdata = mxGetPr(plhs[5]);
        
        plhs[6] = mxCreateNumericArray(3, dims_xS, mxDOUBLE_CLASS, mxREAL);
        dxSdata = mxGetPr(plhs[6]);
        
        dims_yS[0] = nt;
        dims_yS[1] = ny;
        dims_yS[2] = np;
        plhs[7] = mxCreateNumericArray(3, dims_yS, mxDOUBLE_CLASS, mxREAL);
        ySdata = mxGetPr(plhs[7]);
    }
    
    /*
     * ----------------------------------------------------------------
     * IDAS
     * ----------------------------------------------------------------
     */
    
    id = N_VNew_Serial(nx);
    
    idlist = mxGetData(mxGetField(prhs[3], 0, "id"));
    id_tmp = NV_DATA_S(id);
    for (ip=0; ip<np; ip++) {
        id_tmp[ip] = idlist[ip];
    }

    status = IDASetId(ida_mem, id);
    if (status != IDA_SUCCESS) goto error_return;
    
    status = IDASetSuppressAlg(ida_mem, TRUE);
    if (status != IDA_SUCCESS) goto error_return;
    
    /*status = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, RCONST(ts[1]));
    if (status != IDA_SUCCESS) goto error_return;*/
    
    
    ida_status = 0.0;
    
    for (it=0; it < nt; it++) {
        
        /* only integrate if no errors occured */
        if(ida_status == 0.0) {
            /* only integrate after tstart */
            if(ts[it] > tstart) {
                if(nx>0) {
                    ida_status = IDASolve(ida_mem, RCONST(ts[it]), &t, x, dx, IDA_NORMAL);
                }
            }
            tdata[it] = ts[it];
            fu(data, ts[it]);
            fv(data, ts[it], x, dx);
            
            x_tmp = NV_DATA_S(x);
            dx_tmp = NV_DATA_S(dx);
            
            for(ix=0; ix < nx; ix++) {
                xdata[ix*nt+it] = x_tmp[ix];
                dxdata[ix*nt+it] = dx_tmp[ix];
                /* set negative values to zeros */
                if(data->qpositivex[ix]>0.5 && xdata[ix*nt+it]<0.0) xdata[ix*nt+it] = 0.0;
                
            }
        } else {
            for(ix=0; ix < nx; ix++) {
                xdata[ix*nt+it] = 0.0;
                dxdata[ix*nt+it] = 0.0;
            }
        }
        fy(ts[it],nt,it,ydata,data->p,data->k,data->u,xdata);
        
        /* only set output sensitivities if no errors occured */
        if(ida_status == 0.0) {
            if (sensi == 1) {
                if(ts[it] > tstart) {
                    if(nx>0) {
                        status = IDAGetSens(ida_mem, &t, sx);
                        if (status != IDA_SUCCESS) goto error_return;
                        status = IDAGetSensDky(ida_mem, t, 1, sdx);
                        if (status != IDA_SUCCESS) goto error_return;
                    }
                }
                fsu(data, ts[it]);
                fsv(data, ts[it], x, dx);
                
                for(ip=0; ip < np; ip++) {
                    if(nx>0) {
                        sx_tmp = NV_DATA_S(sx[ip]);
                        sdx_tmp = NV_DATA_S(sdx[ip]);
                        for(ix=0; ix < nx; ix++) {
                            xSdata[(ip*nx + ix)*nt + it] = sx_tmp[ix];
                            dxSdata[(ip*nx + ix)*nt + it] = sdx_tmp[ix];
                        }

                    }
                }
                fsy(ts[it],nt,it,plist,nx,ny,np,ySdata,data->p,data->k,data->u,xdata,data->su,xSdata);
            }
        } else {
            if (sensi == 1) {
                for(ip=0; ip < np; ip++) {
                    if(nx>0) {
                        for(ix=0; ix < nx; ix++) {
                            xSdata[(ip*nx + ix)*nt + it] = 0.0;
                        }
                    }
                }
                fsy(ts[it],nt,it,plist,nx,ny,np,ySdata,data->p,data->k,data->u,xdata,data->su,xSdata);
            }
        }
        
    }
    
    
    /*
     * ----------------------------------------------------------------
     * OUTPUT
     * ----------------------------------------------------------------
     */
    
    /*
     *  MATLAB usage: [status,t,x,y,sx,sy,X,S] = simu_ode(tout,theta)
     */
    
    plhs[0] = mxCreateDoubleScalar((double)ida_status);
    
    /* Free memory */
    if(nx>0) {
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(dx);
        if (sensi == 1) {
            N_VDestroyVectorArray_Serial(sx, np);
            N_VDestroyVectorArray_Serial(sdx, np);
        }
        IDAFree(&ida_mem);
    }
    
    N_VDestroy_Serial(id);
    
    mxFree(data);
    mxFree(plist);
    mxFree(p);
    status = ida_status;
    
    return;
    
error_return:
    status = -1;
    plhs[0] = mxCreateDoubleScalar((double)status);
    plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(0,0,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(0,0,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(0,0,mxREAL);
    return;
}