%  This script simulates the PapI regulation model using
%  several modeling approaches and visualizes the simulation results.
%%
clear
clc
close all
modelDefName = 'modelDef_PapI_regulation';
options_plots.save = true;
%% Parameter vector
theta = [1;2.5;-2.25;1;1.2;-0.2;0.01;1.2;-0.2;0.01;2.5;-2.25;10;1]; % from Munsky and Khammash (2006) Table 1 with r/(r+1) ~ 5/6
% Time vector
t = linspace(0,10,100);
%% Simulation using LNA
modelName = 'PapI_LNA';
method = 'LNA';
System_LNA = genSimFile(modelName,modelDefName,method);
System_LNA.sol = simulate_PapI_LNA(t,theta);
plotSSE(System_LNA,options_plots)
%% Simulation using third-order conditional moment equations (MCM3) with zero cumulants closure
modelName = 'PapI_MCM3';
method = 'CMEC_3_ZC_3_a';
System_MCM3 = genSimFileIDA(modelName,modelDefName,method);
System_MCM3.sol = simulate_PapI_MCM3(t,theta);
plotMCM(System_MCM3,options_plots)
%% Simulation using second-order conditional moment equations (MCM2) with low dispersion closure
modelName = 'PapI_MCM2';
method = 'CMEC_2_LD_2_a';
System_MCM2 = genSimFileIDA(modelName,modelDefName,method);
System_MCM2.sol = simulate_PapI_MCM2(t,theta);
plotMCM(System_MCM2,options_plots)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'PapI_MM2';
method = 'MEC_2_LD_2_a';
System_MM2 = genSimFile(modelName,modelDefName,method);
System_MM2.sol = simulate_PapI_MM2(t,theta);
plotMM(System_MM2,options_plots)