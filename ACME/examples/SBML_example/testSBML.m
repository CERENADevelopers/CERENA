clear
clc
close all
importSBML('GeneExpressionCopasi2')
modelName = 'GeneExpressionCopasi2';
modelDefName = 'modelDef_GeneExpressionCopasi2';
%%
t = linspace(0,100,500);
% theta = [0.3;0.3;10;1;4;1;0.015;1;1;0;0;20];
theta = [0.3;0.3;10;1;4;1;0.015];
%% check compilation
method = 'RRE';
System_RRE = genSimFile(modelName,modelDefName,method);

