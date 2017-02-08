%  This script imports the three-stage gene expression model specified in
%  SBML, and simulates the model using reaction rate equations and
%  visualizes the resutls
%%
clear
clc
close all
importSBML('GeneExpressionCopasi')
modelName = 'GeneExpressionCopasi';
modelDefName = 'modelDef_GeneExpressionCopasi';
%%
t = linspace(0,100,500);
% theta = [0.3;0.3;10;1;4;1;0.015;1;1;0;0;20];
theta = [0.3;0.3;10;1;4;1;0.015];
%% check compilation
method = 'RRE';
System_RRE = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
System_RRE.sol = simulate_GeneExpressionCopasi(t,theta);
plotSSE(System_RRE)
