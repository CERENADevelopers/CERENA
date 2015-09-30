clear
close all
clc
%% COMPILATION

[exdir,~,~]=fileparts(which('example_model_2.m'));
cd(exdir)
% compile the model
llhwrap('model_example_2','example_model_2_syms',exdir)

%% SIMULATION

% time vector
t = linspace(0,4,21);
p = [1;0.3;2;3];
k = [];

D.Y = [ 0.0054
    0.0183
   -0.0226
    0.2785
    0.3723
    0.3704
    0.3543
    0.3215
    0.3097
    0.2593
    0.1802
    0.1911
    0.1401
    0.1088
    0.0971
    0.0718
    0.0593
    0.0646
    0.0548
    0.0475
    0.0340];

D.Sigma_Y = 0.01*ones(size(D.Y));


options.sensi = 1;
options.cvode_maxsteps = 1e6;
% load mex into memory
sol = llh_model_example_2(t,log10(p),k,D,options);



figure
errorbar(t,D.Y,D.Sigma_Y)
hold on
plot(t,sol.y)
legend('data','simulation')
xlabel('time t')
ylabel('observable')
title(['log-likelihood: ' num2str(sol.llh) ])

%% FD

eps = 1e-5;
xi = log10(p);
grad_fd = NaN(4,1);
for ip = 1:4;
    options.sensi = 0;
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solp = llh_model_example_2(t,xip,k,D,options);
    grad_fd_f(ip,1) = (solp.llh-sol.llh)/eps;
    xip = xi;
    xip(ip) = xip(ip) - eps;
    solp = llh_model_example_2(t,xip,k,D,options);
    grad_fd_b(ip,1) = -(solp.llh-sol.llh)/eps;
end

figure
scatter(abs(grad_fd_f),abs(sol.sllh))
hold on
scatter(abs(grad_fd_b),abs(sol.sllh))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
plot([1e1,1e4],[1e1,1e4],'k:')
xlim([1e1,1e4])
ylim([1e1,1e4])
legend('forward FD','backward FD')






