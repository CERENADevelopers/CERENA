clear
clc
close all
modelDefName = 'modelDef_geneExpression';
t = linspace(0,100,500);
theta = [0.3;0.3;10;1;4;1;0.015;1;0;1e-10];
kappa = 1;
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'genExp_MM2';
method = 'MEC_2_LD_2_c_f';
System_MM2 = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%%
System_MM2.sol = simulate_AMICI_genExp_timeDep(t,theta,kappa);
options.plot_xo = 1;
options.state_order = 2;
plotMM(System_MM2,options)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'genExp_MCM1';
method = 'CMEC_1_LD_1_c';
System_MCM1 = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
% System_MCM2 = genSimFileIDA(modelName,modelDefName,method);
%%
System_MCM1.sol = simulate_AMICI_genExp_MCM1(t,theta,kappa);
options.plot_xo = 1;
options.state_order = 2;
plotMCM(System_MCM1,options)
