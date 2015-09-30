%  This script simulates the three-stage gene expression model using
%  several modeling approaches and visualizes the simulation results.
%
clear
clc
close all
modelDefName = 'modelDef_geneExpression';
t = linspace(0,30,100);
theta = [0.3;0.3;10;1;4;1;0.015;1;0;0];
kappa = 1;
%% Simulation using RRE
modelName = 'geneExp_RRE';
method = 'RRE';
System_RRE = genSimFile(modelName,modelDefName,method);
System_RRE.sol = simulate_geneExp_RRE(t,theta,kappa);
plotSSE(System_RRE)
%% Simulation using RRE with likelihood calculation
modelName = 'geneExp_RRE_LLH';
method = 'RRE';
t_LLH = 1:10;
data.Y = [1;5;10;14;18;21;23;25;26;27];
data.Sigma_Y = 0.1*ones(length(t_LLH),1);
System_RRE = genSimFileLLH(modelName,modelDefName,method);
System_RRE.sol = llh_geneExp_RRE_LLH(t_LLH,theta,kappa,data);
%% Simulation using LNA
modelName = 'geneExp_LNA';
method = 'LNA';
System_LNA = genSimFile(modelName,modelDefName,method);
System_LNA.sol = simulate_geneExp_LNA(t,theta,kappa);
plotSSE(System_LNA)
%% Simulation using EMRE
modelName = 'geneExp_EMRE';
method = 'EMRE';
System_EMRE = genSimFile(modelName,modelDefName,method);
System_EMRE.sol = simulate_geneExp_EMRE(t,theta,kappa);
plotSSE(System_EMRE)
%% Simulation using IOS
modelName = 'geneExp_IOS';
method = 'IOS';
System_IOS = genSimFile(modelName,modelDefName,method);
System_IOS.sol = simulate_geneExp_IOS(t,theta,kappa);
options.lw = 2;
options.fs = 16;
plotSSE(System_IOS,options)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'geneExp_MM2';
method = 'MEC_2_LD_2_c';
System_MM2 = genSimFile(modelName,modelDefName,method);
System_MM2.sol = simulate_geneExp_MM2(t,theta,kappa);
options.plot_xo = 1;
options.state_order = 2;
plotMM(System_MM2,options)
%% Simulation using third-order conditional moment equations (MCM3) with zero cumulants closure
modelName = 'geneExp_MCM2';
method = 'CMEC_3_ZC_2_a';
System_MCM3 = genSimFileIDA(modelName,modelDefName,method);
System_MCM3.sol = simulate_geneExp_MCM2(t,theta,kappa);
plotMCM(System_MCM3)
%% Providing symbolic and numeric initial conditions before compilation of MEX files
modelName = 'geneExp_RRE_IC';
method = 'RRE';
System_RRE_IC = genmexp(modelName,modelDefName,method);
syms r0 p0
ICsym = [1;0;r0;p0;zeros(10,1)];
ICnum = [1;0;0;20;zeros(10,1)];
IC = ICnum;
cvodewrap(modelName,[method,'_',modelName,'_syms'],[],IC)
System_RRE_IC.sol = simulate_geneExp_RRE_IC(t,theta,kappa);
plotSSE(System_RRE_IC)
IC = ICsym;
cvodewrap(modelName,[method,'_',modelName,'_syms'],[],IC)
System_RRE_IC.sol = simulate_geneExp_RRE_IC(t,theta,kappa);
plotSSE(System_RRE_IC)
%% Numeric simulation with arbitrary initial conditions, i.e. different from the IC provided in modeDef file
modelName = 'geneExp_RRE';
method = 'RRE';
System_RRE = genSimFile(modelName,modelDefName,method);
IC = [1;0;0;20;zeros(10,1)];
ind_IC = ones(size(IC));
kappa_IC = [kappa;ind_IC;IC];
System_RRE.sol = simulate_geneExp_RRE(t,theta,kappa_IC);
plotSSE(System_RRE)
%% Simulation using FSP
System_FSP = simulateFSP_matlab(modelDefName,t,theta,kappa);
plotFSP_CERENA(System_FSP)
%% Simulation using SSA
Nssa = 20;
options.mode = 'constant';
options.scale = 'concentration';
System_SSA = simulateSSA_matlab(modelDefName,t,theta,kappa,Nssa,options);
plotSSA(System_SSA)
%% Calculation and visualization of correlation and partial correlation maps
% Correlation
options.visualization = 'on';
options.movieName = 'corr_geneExp';
[corrMat,corrMatAll,covMat] = corrmat(System_MM2,System_MM2.sol.x,options);
% Partial correlation
options.movieName = 'pcorr_geneExp';
[pcorrMat,pcorrMatAll] = pcorrmat(System_MM2,System_MM2.sol.x,covMat,corrMatAll,options);

