clear 
clc
close all
modelDefName = 'modelDef_geneCascade';
%% Parameter vector
theta = [50;0.8;0.04;15;1];
t = linspace(0,0.1,100);
%% Simulation using FSP
% Run the FSP simulation - takes time!
eval(modelDefName)
System.state.mu0 = subs(System.state.mu0,System.parameter.variable,theta);
System_FSP = completeSystemFSP(System);
System_FSP = simulate_FSP(System_FSP,t,theta);
plotFSP_ACME(System_FSP)
% Or load the FSP simulation results
load('FSP_results.mat')
%% Central moment equations, absolute scale, 2nd order, with derivative matching closure
modelName_MM2 = 'geneCascade_MM2';
method = 'MEC_2_DM_1_a';
System_MM2 = genSimFile(modelName_MM2,modelDefName,method);
System_MM2.sol = simulate_geneCascade_MM2(t,theta);
plotMM(System_MM2)