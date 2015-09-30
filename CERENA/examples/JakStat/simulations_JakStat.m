%  This script simulates the Jak-Stat signalling pathway using
%  several modeling approaches and visualizes the simulation results.
%%
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
options_plots.save = false;
%% Simulation using RRE
modelName = 'JakStat_RRE';
method = 'RRE';
% System_RRE = genSimFile(modelName,modelDefName,method);
System_RRE = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('RRE_JakStat_RRE_syms.m');
clear('RRE_JakStat_RRE_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('RRE_JakStat_RRE_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_RRE = load('JakStat_RRE_RRE_System');
System_RRE.sol = simulate_JakStat_RRE(t,theta,kappa);
plotSSE(System_RRE,options_plots)
%% Simulation using LNA
modelName = 'JakStat_LNA';
method = 'LNA';
% System_LNA = genSimFile(modelName,modelDefName,method);
System_LNA = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('LNA_JakStat_LNA_syms.m');
clear('LNA_JakStat_LNA_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('LNA_JakStat_LNA_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_LNA = load('JakStat_LNA_LNA_System');
System_LNA.sol = simulate_JakStat_LNA(t,theta,kappa);
plotSSE(System_LNA,options_plots)
%% Simulation using EMRE
modelName = 'JakStat_EMRE';
method = 'EMRE';
% System_EMRE = genSimFile(modelName,modelDefName,method);
System_EMRE = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('EMRE_JakStat_EMRE_syms.m');
clear('EMRE_JakStat_EMRE_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('EMRE_JakStat_EMRE_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_EMRE = load('JakStat_EMRE_EMRE_System');
System_EMRE.sol = simulate_JakStat_EMRE(t,theta,kappa);
plotSSE(System_EMRE,options_plots)
%% Simulation using IOS
modelName = 'JakStat_IOS';
method = 'IOS';
% System_IOS = genSimFile(modelName,modelDefName,method);
System_IOS = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('IOS_JakStat_IOS_syms.m');
clear('IOS_JakStat_IOS_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('IOS_JakStat_IOS_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_IOS = load('JakStat_IOS_IOS_System');
System_IOS.sol = simulate_JakStat_IOS(t,theta,kappa);
options_plots.plot_xsse = false;
plotSSE(System_IOS,options_plots)
%% Simulation using second-order moment equations (MM2) with low dispersion closure
modelName = 'JakStat_MM2';
method = 'MEC_2_LD_2_c';
% System_MM2 = genSimFile(modelName,modelDefName,method);
System_MM2 = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('MEC_2_LD_2_c_JakStat_MM2_syms.m');
clear('MEC_2_LD_2_c_JakStat_MM2_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('MEC_2_LD_2_c_JakStat_MM2_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_MM2 = load('JakStat_MM2_MEC_2_LD_2_c_System');
System_MM2.sol = simulate_JakStat_MM2(t,theta,kappa);
plotMM(System_MM2,options_plots)
%% Simulation using third-order moment equations (MM3) with low dispersion closure
modelName = 'JakStat_MM3';
method = 'MEC_3_LD_3_c';
% System_MM3 = genSimFile(modelName,modelDefName,method);
System_MM2 = genmexp(modelName,modelDefName,method);
% Replace the input by the spline function
s = fileread('MEC_3_LD_3_c_JakStat_MM3_syms.m');
clear('MEC_3_LD_3_c_JakStat_MM3_syms.m');
s = strrep(s,'pEpoR','spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)');
fid = fopen('MEC_3_LD_3_c_JakStat_MM3_syms.m','w');
fprintf(fid,'%s',s);
fclose(fid);
rehash
cvodewrap(modelName,[method,'_',modelName,'_syms'])
% System_MM3 = load('JakStat_MM3_MEC_3_LD_3_c_System');
System_MM3.sol = simulate_JakStat_MM3(t,theta,kappa);
plotMM(System_MM3,options_plots)
%% Simulation using SSA %% SSA
% theta_SSA = [p;1;840;270;scalePar;sp;0.2202;0.0746];
eval(modelDefName);
System.scaleConversion = 'Macro_to_Micro';
System_SSA = completeSystem(System);
System_SSA = completeSystemSSA(System_SSA);
System_SSA.v = @(t,x,theta)[spline_pos5_matlab(t, 0.0, theta(12), 5.0, theta(13), 10.0, theta(14), 20.0, theta(15), 60.0, theta(16), 0, 0.0)*theta(1)*x(1);...
    -(theta(2)*(x(2)/2-x(2)^2/2))/(theta(5)*theta(6));...
    theta(3)*x(3);...
    theta(4)*x(4);...
    theta(4)*x(5);...
    theta(4)*x(6);...
    theta(4)*x(7);...
    theta(4)*x(8);...
    theta(4)*x(9)];
options.mode = 'time-dependent';
options.scale = 'concentration';
Nssa = 1;
tic
System_SSA.sol = simulate_SSA(System_SSA,t,theta,kappa,Nssa,options);
toc
System_SSA.sol.theta = theta;
System_SSA.sol.kappa = kappa;
plotSSA(System_SSA,options_plots)
% Or load the results of 1000 SSA runs
% load('1000SSA_SensParam_spline')
%% Calculation and visualization of correlation and partial correlation maps
% Correlation
options.visualization = 'on';
options.movieName = 'corr_JakStat';
[corrMat,corrMatAll,covMat] = corrmat(System_MM2,System_MM2.sol.x,options);
% Partial correlation
options.movieName = 'pcorr_JakStat';
[pcorrMat,pcorrMatAll] = pcorrmat(System_MM2,System_MM2.sol.x,covMat,corrMatAll,options);