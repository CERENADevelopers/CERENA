#include "geneExp_LNA.c"
                
                int cvodewrap_init(void *cvode_mem, N_Vector x, double t){
                    return CVodeInit(cvode_mem, xdot_geneExp_LNA, RCONST(t), x);
                }
                int cvodewrap_binit(void *cvode_mem, int which, N_Vector xB, double t){
                    return CVodeInitB(cvode_mem, which, xBdot_geneExp_LNA, RCONST(t), xB);
                }
                int cvodewrap_qbinit(void *cvode_mem, int which, N_Vector xQB){
                    return CVodeQuadInitB(cvode_mem, which, xQB_geneExp_LNA, xQB);
                }
                
                void fx0(N_Vector x0, void *user_data){
                    UserData data = (UserData) user_data;
                    x0_geneExp_LNA(x0, data);
                }
                
                int cvodewrap_SetDenseJacFn(void *cvode_mem){
                    return CVDlsSetDenseJacFn(cvode_mem, J_geneExp_LNA);
                }
                int cvodewrap_SetSparseJacFn(void *cvode_mem){
                    return CVSlsSetSparseJacFn(cvode_mem, JSparse_geneExp_LNA);
                }
                int cvodewrap_SetBandJacFn(void *cvode_mem){
                    return CVDlsSetBandJacFn(cvode_mem, JBand_geneExp_LNA);
                }
                int cvodewrap_SetJacTimesVecFn(void *cvode_mem){
                    return CVSpilsSetJacTimesVecFn(cvode_mem, Jv_geneExp_LNA);
                }
                int cvodewrap_SetDenseJacFnB(void *cvode_mem,int which){
                    return CVDlsSetDenseJacFnB(cvode_mem, which, JB_geneExp_LNA);
                }
                int cvodewrap_SetSparseJacFnB(void *cvode_mem, int which){
                    return CVSlsSetSparseJacFnB(cvode_mem, which, JSparseB_geneExp_LNA);
                }
                int cvodewrap_SetBandJacFnB(void *cvode_mem,int which){
                    return CVDlsSetBandJacFnB(cvode_mem, which, JBBand_geneExp_LNA);
                }
                int cvodewrap_SetJacTimesVecFnB(void *cvode_mem,int which){
                    return CVSpilsSetJacTimesVecFnB(cvode_mem, which, JvB_geneExp_LNA);
                }
                int fJ(long int N, realtype t, N_Vector x,N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
                    return J_geneExp_LNA(N,t,x,fx,J,user_data,tmp1,tmp2,tmp3);
                }
                int fJB(long int N, realtype t, N_Vector x, N_Vector xB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
                    return JB_geneExp_LNA(N,t,x,xB,fx,J,user_data,tmp1,tmp2,tmp3);
                }
                
                void fsx0(N_Vector *sx0, void *user_data){
                    UserData data = (UserData) user_data;
                    int ip;
                    int *plist = data->plist;
                    int np = *data->np;
                    for (ip=0; ip<np; ip++) {
                       sx0_geneExp_LNA(plist[ip], sx0[ip], data);
                    }
                }
                
                int cvodewrap_SensInit1(void *cvode_mem, int np, int sensi_meth, N_Vector *sx){
                    return CVodeSensInit1(cvode_mem, np, sensi_meth, sx_geneExp_LNA, sx);
                }
                int cvodewrap_RootInit(void *cvode_mem, int nr){
                    return CVodeRootInit(cvode_mem, nr, root_geneExp_LNA);
                }
                void froot(double t, N_Vector x, realtype *gout, void *user_data){
                    root_geneExp_LNA(t, x, gout, user_data);
                }
                
                void fy(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
                    y_geneExp_LNA(t, nt, it, y, p, k, u, x);
                }
                
                void fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data){
                    xdot_geneExp_LNA(t,x,xdot,user_data);
                }
                void fdydp(double t, int nt, int it,double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
                    dydp_geneExp_LNA(t, nt, it, dydp, y, p, k, u, x, plist, np, ny);
                }
                void fdydx(double t,double *dydx, double *y, double *p, double *k, double *x){
                    dydx_geneExp_LNA(t, dydx, y, p, k, x);
                }
                
                
                void fsy(double t, int nt, int it, int *plist, int nx, int ny, int np, double *sy, double *p, double *k, double *x, double *sx){
                    int ip;
                    for (ip=0; ip<np; ip++) {
                        sy_geneExp_LNA(t, nt, it, plist[ip], ip,  nx, ny,  sy, p, k, x, sx);
                    }
                }
                double fsroot(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(sroot_geneExp_LNA(t, ip, ir, x, sx, user_data));
                }
                double fsrootval(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(srootval_geneExp_LNA(t, ip, ir, x, sx, user_data));
                }
                double fs2root(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(s2root_geneExp_LNA(t, ip, jp, ir, x, sx, user_data));
                }
                double fs2rootval(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
                    return(s2rootval_geneExp_LNA(t, ip, jp, ir, x, sx, user_data));
                }
                void cvodewrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data){
                    char buffer [250];
                    sprintf(buffer,"CVODES ERROR: in module %s in function %s : %s ",module,function,msg);
                    mexWarnMsgTxt(buffer);
                }

