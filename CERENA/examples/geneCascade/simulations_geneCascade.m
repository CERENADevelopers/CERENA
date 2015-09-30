%  This script simulates the gene cascade model using
%  several modeling approaches and visualizes the simulation results.
%%
clear 
clc
close all
modelDefName = 'modelDef_geneCascade';
%% Parameter vector
theta = [50;0.8;0.04;15;1];
t = linspace(0,0.1,100);
%% Simulation using FSP
% Run the FSP simulation - this takes a while!
System_FSP = simulateFSP_matlab(modelDefName,t,theta);
plotFSP_CERENA(System_FSP)
% Or load the FSP simulation results
% load('FSP_results.mat')
%% Central moment equations, absolute scale, 2nd order, with derivative matching closure
modelName_MM2 = 'geneCascade_MM2';
method = 'MEC_2_DM_1_a';
System_MM2 = genSimFile(modelName_MM2,modelDefName,method);
System_MM2.sol = simulate_geneCascade_MM2(t,theta);
plotMM(System_MM2)