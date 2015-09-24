#include "geneExp_MCM3.c"
                
                int idawrap_init(void *ida_mem, N_Vector x, N_Vector dx, double t){
                    return IDAInit(ida_mem, fx_geneExp_MCM3, RCONST(t), x, dx);
                }
                int idawrap_binit(void *ida_mem, int which, N_Vector xB, N_Vector dxB, double t){
                    return IDAInitB(ida_mem, which, fxB_geneExp_MCM3, RCONST(t), xB, dxB);
                }
                int idadewrap_qbinit(void *ida_mem, int which, N_Vector xQB){
                    return IDAQuadInitB(ida_mem, which, fxQB_geneExp_MCM3, xQB);
                }
                
                void fx0(N_Vector x0, void *user_data){
                    UserData data = (UserData) user_data;
                    fx0_geneExp_MCM3(x0, data);
                }
                
                void fdx0(N_Vector x0, N_Vector dx0, void *user_data){
                    UserData data = (UserData) user_data;
                    fdx0_geneExp_MCM3(x0, dx0, data);
                }
                
                int idawrap_SetDenseJacFn(void *ida_mem){
                    return IDADlsSetDenseJacFn(ida_mem, DenseJacFn_geneExp_MCM3);
                }
                
                void fsx0(int is, N_Vector sx_is, void *user_data){
                    UserData data = (UserData) user_data;
                    fsx0_geneExp_MCM3(is, sx_is, data);
                }
                void fsdx0(int is, N_Vector sdx_is, void *user_data){
                    UserData data = (UserData) user_data;
                    fsdx0_geneExp_MCM3(is, sdx_is, data);
                }
                
                int idawrap_SensInit(void *ida_mem, int nps, int sensi_meth, N_Vector *sx, N_Vector *sdx){
                    return IDASensInit(ida_mem, nps, sensi_meth, fsx_geneExp_MCM3, sx, sdx);
                }
                
                void fu(void *user_data, double t){
                    UserData data = (UserData) user_data;
                    fu_geneExp_MCM3(data, t);
                }
                
                void fsu(void *user_data, double t){
                    UserData data = (UserData) user_data;
                    fsu_geneExp_MCM3(data, t);
                }
                
                void fv(void *user_data, double t, N_Vector x, N_Vector dx){
                    UserData data = (UserData) user_data;
                    fv_geneExp_MCM3(t, x, dx, data);
                }
                
                void fsv(void *user_data, double t, N_Vector x, N_Vector dx){
                    UserData data = (UserData) user_data;
                    dvdu_geneExp_MCM3(t, x, dx, data);
                    dvdx_geneExp_MCM3(t, x, dx, data);
                    
                }
                
                void fy(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
                    fy_geneExp_MCM3(t, nt, it, y, p, k, u, x);
                }
                void dfydp(double t, int nt, int it,double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int nps, int nobs){
                    dfydp_geneExp_MCM3(t, nt, it, dydp, y, p, k, u, x, plist, nps, nobs);
                }
                void dfydx(double t,double *dydx, double *y, double *p, double *k, double *u, double *x){
                    dfydx_geneExp_MCM3(t, dydx, y, p, k, u, x);
                }
                
                void fsy(double t, int nt, int it, int *plist, int nx, int ny, int nps, double *sy, double *p, double *k, double *u, double *x, double *su, double *sx){
                    int ip;
                    for (ip=0; ip<nps; ip++) {
                        fsy_geneExp_MCM3(t, nt, it, plist[ip], ip,  nx, ny,  sy, p, k, u, x, su, sx);
                    }
                }
