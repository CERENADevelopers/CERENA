clear 
clc
close all
modelDefName = 'modelDef_JakStat';
sp = 10.^[-2.7859; -0.2561;-0.0751;-0.4112;-5.0000];
scalePar = 10.^[-0.6415;-0.7353;0.0273;-0.1079];
p = 10.^[0.5951;1;-0.9489;-0.0075];
theta = [p;1;840;270;scalePar;sp];
kappa = [];
t = linspace(0,60,200);
%% Simulation using RRE
modelName = 'JakStat_RRE';
method = 'RRE';
System_RRE = genSimFile(modelName,modelDefName,method);
System_RRE.sol = simulate_JakStat_RRE(t,theta,kappa);
plotSSE(System_RRE)
%% Simulation using LNA
modelName = 'JakStat_LNA';
method = 'LNA';
System_LNA = genSimFile(modelName,modelDefName,method);
System_LNA.sol = simulate_JakStat_LNA(t,theta,kappa);
plotSSE(System_LNA)
%% Simulation using EMRE
modelName = 'JakStat_EMRE';
method = 'EMRE';
System_EMRE = genSimFile(modelName,modelDefName,method);
System_EMRE.sol = simulate_JakStat_EMRE(t,theta,kappa);
plotSSE(System_EMRE)
%% Simulation using IOS
modelName = 'JakStat_IOS';
method = 'IOS';
System_IOS = genSimFile(modelName,modelDefName,method);
System_IOS.sol = simulate_JakStat_IOS(t,theta,kappa);
plotSSE(System_IOS)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'JakStat_MM2';
method = 'MEC_2_LD_2_c';
System_MM2 = genSimFile(modelName,modelDefName,method);
System_MM2.sol = simulate_JakStat_MM2(t,theta,kappa);
plotMM(System_MM2)
%% Simulation using third-order moment equations (MM3) with low dispersion closure
modelName = 'JakStat_MM3';
method = 'MEC_3_LD_3_c';
System_MM3 = genSimFile(modelName,modelDefName,method);
System_MM3.sol = simulate_JakStat_MM3(t,theta,kappa);
plotMM(System_MM3)
%% Simulation using SSA %% SSA
theta_SSA = [p;1;840;270;scalePar;sp;0.2202;0.0746];
eval(modelDefName);
System_SSA = completeSystem(System);
System_SSA = completeSystemSSA(System_SSA);
options.mode = 'time-dependent';
options.scale = 'concentration';
Nssa = 2;
tic
System_SSA.sol = simulate_SSA(System_SSA,t,theta_SSA,kappa,Nssa,options);
System_SSA.sol.theta = theta;
System_SSA.sol.kappa = kappa;
toc
plotSSA(System_SSA)
% Or load the results of 1000 SSA runs
load('1000SSA_SensParam')
%% Calculation and visualization of correlation and partial correlation maps
% Correlation
options.visualization = 'on';
options.movieName = 'corr_JakStat';
[corrMat,corrMatAll,covMat] = corrmat(System_MM2,System_MM2.sol.x,options);
% Partial correlation
options.movieName = 'pcorr_JakStat';
[pcorrMat,pcorrMatAll] = pcorrmat(System_MM2,System_MM2.sol.x,covMat,corrMatAll,options);