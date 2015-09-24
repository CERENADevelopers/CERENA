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
plotSSE(System_IOS)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'geneExp_MM2';
method = 'MEC_2_LD_2_c';
System_MM2 = genSimFile(modelName,modelDefName,method);
System_MM2.sol = simulate_geneExp_MM2(t,theta,kappa);
plotMM(System_MM2)
%% Simulation using third-order conditional moment equations (MCM3) with zero cumulants closure
modelName = 'geneExp_MCM3';
method = 'CMEC_3_ZC_3_a';
System_MCM3 = genSimFileIDA(modelName,modelDefName,method);
System_MCM3.sol = simulate_geneExp_MCM3(t,theta,kappa);
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
eval(modelDefName)
System.state.mu0 = subs(System.state.mu0,System.parameter.variable,theta);
System.state.mu0 = subs(System.state.mu0,System.kappa.variable,kappa);
System_FSP = completeSystemFSP(System);
theta_FSP = [theta;kappa];
System_FSP = simulate_FSP(System_FSP,t,theta_FSP);
plotFSP_ACME(System_FSP)
%% Simulation using SSA
eval(modelDefName);
System_SSA = completeSystem(System);
System_SSA = completeSystemSSA(System_SSA);
options.mode = 'constant';
options.scale = 'concentration';
Nssa = 20;
tic
System_SSA.sol = simulate_SSA(System_SSA,t,theta,kappa,Nssa,options);
System_SSA.sol.theta = theta;
System_SSA.sol.kappa = kappa;
toc
plotSSA(System_SSA)
%% Calculation and visualization of correlation and partial correlation maps
% Correlation
options.visualization = 'on';
options.movieName = 'corr_geneExp';
[corrMat,corrMatAll,covMat] = corrmat(System_MM2,System_MM2.sol.x,options);
% Partial correlation
options.movieName = 'pcorr_geneExp';
[pcorrMat,pcorrMatAll] = pcorrmat(System_MM2,System_MM2.sol.x,covMat,corrMatAll,options);

