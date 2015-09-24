clear
clc
close all
modelDefName = 'modelDef_PapI_regulation';
%% Parameter vector
theta = [1;2.5;-2.25;1;1.2;-0.2;0.01;1.2;-0.2;0.01;2.5;-2.25;10;1]; % from Munsky and Khammash (2006) Table 1 with r/(r+1) ~ 5/6
% Time vector
t = linspace(0,10.5,100);
%% Simulation using LNA
modelName = 'PapI_LNA';
method = 'LNA';
System_LNA = genSimFile(modelName,modelDefName,method);
System_LNA.sol = simulate_PapI_LNA(t,theta);
plotSSE(System_LNA)
%% Simulation using third-order conditional moment equations (MCM3) with low dispersion closure
modelName = 'PapI_MCM3';
method = 'CMEC_3_ZC_3_a';
System_MCM3 = genSimFileIDA(modelName,modelDefName,method);
System_MCM3.sol = simulate_PapI_MCM3(t,theta);
plotMCM(System_MCM3)
