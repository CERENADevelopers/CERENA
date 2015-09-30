%  This script simulates the Jak-Stat signaling pathway using
%  the second-order moment equations and calculates the negative
%  log-likelihood values for given data. In addition, it calculates the
%  forward and adjoint sensitivities.
%%
clear
clc
close all
modelDefName = 'modelDef_JakStat';
%% Parameter vector
% With spline input
sp = 10.^[-2.7859; -0.2561;-0.0751;-0.4112;-5.0000];
scalePar = 10.^[-0.6415;-0.7353;0.0273;-0.1079];
p = 10.^[0.5951;1;-0.9489;-0.0075];
theta = [p;1;840;270;scalePar;sp];
kappa = [];
% Time vector
t = [0;2;4;6;8;10;12;14;16;18;20;25;30;40;50;60];
% Data to be used in llhwrap to calculate the objective function value
data.Y = [0,1;0.3315,0.9275;0.8645,0.7923;0.9635,0.7778;0.9279,0.7053;0.8162,0.6522;...
          0.7553,0.5894;0.7680,0.5894;0.8416,0.6377;0.7680,0.6425;0.8010,0.6908;0.7832,0.6908;...
          0.8086,0.7585;0.4888,0.8068;0.2782,0.9275;0.2553,0.9710];
data.Sigma_Y = ones(length(t),2);
%% Central moment equations, concentration scale, 2nd order, compiled with cvodewrap
modelName_cw = 'JakStat_MM2_cw';
method = 'MEC_2_LD_1_c';
% System_MM2_cw = genSimFile(modelName_cw,modelDefName,method);
System_MM2_cw = genmexp(modelName_cw,modelDefName,method);
% Replace the input by the spline function
s = fileread('MEC_2_LD_1_c_JakStat_MM2_cw_syms.m');
clear('MEC_2_LD_1_c_JakStat_MM2_cw_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('MEC_2_LD_1_c_JakStat_MM2_cw_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName_cw,[method,'_',modelName_cw,'_syms'])
options.sensi = 1;
options.sensi_meth = 1; % Forward sensitivity analysis
% Forward sensitivities for all state variables
System_MM2_cw.sol = simulate_JakStat_MM2_cw(t,theta,kappa,options);
options.sensi_meth = 2; % Adjoint sensitivity analysis
% Adjoint sensitivities for all output variables
System_MM2_cw.sol = simulate_JakStat_MM2_cw(t,theta,kappa,options);

%% Central moment equations, concentration scale, 2nd order, compiled with llhwrap
modelName_lw = 'JakStat_MM2_lw';
% System_MM2_lw = genSimFileLLH(modelName_lw,modelDefName,method);
System_MM2_lw = genmexp(modelName_lw,modelDefName,method);
% Replace the input by the spline function
s = fileread('MEC_2_LD_1_c_JakStat_MM2_lw_syms.m');
clear('MEC_2_LD_1_c_JakStat_MM2_lw_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('MEC_2_LD_1_c_JakStat_MM2_lw_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
llhwrap(modelName_lw,[method,'_',modelName_lw,'_syms'])
options.sensi_meth = 2; % Adjoint sensitivities
% Adjoint sensitivities for the objective funciton value
System_MM2_lw.sol = llh_JakStat_MM2_lw(t,theta,kappa,data,options);