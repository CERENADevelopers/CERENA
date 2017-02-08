clear
clc
close all
modelDefName = 'modelDef_geneExpression';
t = linspace(0,100,500);
theta = [0.3;0.3;10;1;4;1;0.015;1;0;1e-10];
kappa = 1;
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'geneExpressionIOS';
method = 'IOS';
System_IOS = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
System_IOS.sol = simulate_geneExpressionIOS(t,theta,kappa);
plotSSE(System_IOS)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'geneExpressionMM';
method = 'MEC_3_LD_2_c_f';
System_MM = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%%
System_MM.sol = simulate_geneExpressionMM(t,theta,kappa);
options.plot_xo = 1;
options.state_order = 2;
plotMM(System_MM,options)
%%
modelName = 'geneExpressionsMA'; % The name that will be used for naming the simulation files
method = 'MEC_2_LD_2_c_g'; % Specifying the modeling approach
System_sMA.reduction_order = 1;
System_sMA = genmexp(modelName,modelDefName,method,System_sMA);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%%
System_sMA.sol = simulate_geneExpressionsMA(t,theta,kappa);
options.plot_xo = 1;
options.state_order = 2;
plotMM(System_sMA,options)

%% Simulation using second-order moment equations (MM2) with low dispersion closure
% modelName = 'geneExpressionMCM';
% method = 'CMEC_1_LD_1_c';
% System_MCM = genmexp(modelName,modelDefName,method);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% %%
% [status_MCM,tout_MCM,x_MCM,dx_MCM,mx_MCM,y_MCM] = simulate_geneExpressionMCM(t,theta,kappa);
% %%
% System_MCM.sol = simulate_geneExpressionMCM(t,theta,kappa);
% options.plot_xo = 1;
% options.state_order = 2;
% plotMCM(System_MCM,options)
