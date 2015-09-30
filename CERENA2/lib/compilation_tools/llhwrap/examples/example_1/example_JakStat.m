clear
close all
% compile the model
[exdir,~,~]=fileparts(which('example_JakStat.m'));
llhwrap('model_example_1','JakStat_syms',exdir)

num = xlsread('pnas_data_original.xls');
load('workspace')

t = num(:,1);

D.Y = num(:,[2,4,6]);
D.Sigma_Y = ones(size(t))*(10.^(ar.p(10:12)));


kappa = [1.4,0.45];

xi = ar.p([4 5 6 7 1 13 14 15 16 17 3 2 9 8]);

options.sensi = 0;
sol = llh_model_example_1(t,xi,kappa,D,options);

figure
for iy = 1:3
    subplot(2,2,iy)
    errorbar(t,D.Y(:,iy),D.Sigma_Y(:,iy),'LineWidth',1.5)
    hold on
    plot(t,sol.y(:,iy),'.-','LineWidth',1.5)
    xlim([0,60])
    xlabel('t')
    switch(iy)
        case 1
            ylabel('pStat')
        case 2
            ylabel('tStat')
        case 3
            ylabel('pEpoR')
    end
    ylim([0,1.2])
    set(gca,'FontSize',15)
    set(gca,'LineWidth',1.5)
end

% generate random parameter
xi_rand = xi + 0.1;
options.sensi = 1;
sol = llh_model_example_1(t,xi_rand,kappa,D,options);

options.sensi = 0;
eps = 1e-4;
fd_grad = NaN(length(xi),1);
for ip = 1:length(xi)
    xip = xi_rand;
    xip(ip) = xip(ip) + eps;
    psol = llh_model_example_1(t,xip,kappa,D,options);
    fd_grad(ip) = (psol.llh-sol.llh)/eps;
end

figure
scatter(abs(sol.sllh),abs(fd_grad))
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1e0,1e4])
ylim([1e0,1e4])
box on
set(gca,'FontSize',15)
set(gca,'LineWidth',1.5)
hold on
plot([1e0,1e4],[1e0,1e4],'k:')
ylabel('absolute value FD gradient entries')
xlabel('absolute value AS gradient entries')

options.sensi = 1;
options.ordering = 0;
while true
    sol = llh_model_example_1(t,xi,kappa,D,options);
end

