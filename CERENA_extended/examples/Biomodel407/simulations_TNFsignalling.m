clear
clc
close all
modelDefName = 'modelDef_biomodel407';
load('theta_curated.mat') % load parameter values
load('kappa_curated.mat') % load initial conditions
t = linspace(0,50000,10000);
%% Derivation of MA and sMA equations of degree 1, 2 and 3;
%%%% This will take a while; instead you can use the pregenerated models
%%%% below.
%
% %% full MM simulation
% modelName = 'biomodel407_MM2_f_a';
% method = 'MEC_2_LD_1_a_f_2';
% System_MM2_f = genmexp(modelName,modelDefName,method);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_f done')
% %% reduced MM with zero closure
% modelName = 'biomodel407_MM2_g_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g = genmexp(modelName,modelDefName,method);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g done')
% %% reduced MM order 2 with zero closure
% modelName = 'biomodel407_MM2_g2_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g2.reduction_order = 2;
% System_MM2_g2 = genmexp(modelName,modelDefName,method,System_MM2_g2);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g2 done')
% %% reduced MM order 3 with zero closure
% modelName = 'biomodel407_MM2_g3_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g3.reduction_order = 3;
% System_MM2_g3 = genmexp(modelName,modelDefName,method,System_MM2_g3);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g3 done')

%%
load('biomodel407_MM2_f_a_MEC_2_LD_1_a_f_2_System.mat')
System_MM2_f = System;
load('biomodel407_MM2_g_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g = System;
load('biomodel407_MM2_g2_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g2 = System;
load('biomodel407_MM2_g3_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g3 = System;
%%
System_MM2_f.sol = simulate_biomodel407_MM2_f_a(t,theta,kappa);
System_MM2_g.sol = simulate_biomodel407_MM2_g_a(t,theta,kappa);
System_MM2_g2.sol = simulate_biomodel407_MM2_g2_a(t,theta,kappa);
System_MM2_g3.sol = simulate_biomodel407_MM2_g3_a(t,theta,kappa);

ind_mean_f = find(sum(System_MM2_f.MM.sym.state.order>=1,2) == 1);
ind_mean_g = find(sum(System_MM2_g.MM.sym.state.order>=1,2) == 1);
ind_mean_g2 = find(sum(System_MM2_g2.MM.sym.state.order>=1,2) == 1);
ind_mean_g3 = find(sum(System_MM2_g3.MM.sym.state.order>=1,2) == 1);

ind_var_f = find((System_MM2_f.MM.sym.state.order(:,end-1)== System_MM2_f.MM.sym.state.order(:,end)).*(sum(System_MM2_f.MM.sym.state.order~=0,2)==2));
ind_var_g = find((System_MM2_g.MM.sym.state.order(:,end-1)== System_MM2_g.MM.sym.state.order(:,end)).*(sum(System_MM2_g.MM.sym.state.order~=0,2)==2));
ind_var_g2 = find((System_MM2_g2.MM.sym.state.order(:,end-1)== System_MM2_g2.MM.sym.state.order(:,end)).*(sum(System_MM2_g2.MM.sym.state.order~=0,2)==2));
ind_var_g3 = find((System_MM2_g3.MM.sym.state.order(:,end-1)== System_MM2_g3.MM.sym.state.order(:,end)).*(sum(System_MM2_g3.MM.sym.state.order~=0,2)==2));

mean_x_f = System_MM2_f.sol.x(:,ind_mean_f);
var_x_f = System_MM2_f.sol.x(:,ind_var_f);
mean_x_g = System_MM2_g.sol.x(:,ind_mean_g);
var_x_g = System_MM2_g.sol.x(:,ind_var_g);
mean_x_g2 = System_MM2_g2.sol.x(:,ind_mean_g2);
var_x_g2 = System_MM2_g2.sol.x(:,ind_var_g2);
mean_x_g3 = System_MM2_g3.sol.x(:,ind_mean_g3);
var_x_g3 = System_MM2_g3.sol.x(:,ind_var_g3);

e_mean_x_g = bsxfun(@rdivide,abs(mean_x_g-mean_x_f),max(mean_x_f));
e_var_x_g = bsxfun(@rdivide,abs(var_x_g-var_x_f),max(var_x_f));
e_mean_x_g2 = bsxfun(@rdivide,abs(mean_x_g2-mean_x_f),max(mean_x_f));
e_var_x_g2 = bsxfun(@rdivide,abs(var_x_g2-var_x_f),max(var_x_f));
e_mean_x_g3 = bsxfun(@rdivide,abs(mean_x_g3-mean_x_f),max(mean_x_f));
e_var_x_g3 = bsxfun(@rdivide,abs(var_x_g3-var_x_f),max(var_x_f));

%% Visualization
fs=7;
lw = 1;
ms = 2;

syms pCasp3 pCasp8 TNFR_E NFkB
ind_species_plot= find(ismember(System_MM2_f.state.variable,[pCasp8; TNFR_E; NFkB]));
yt = [0:0.003:0.006;0:0.4:0.8;0:0.12:0.24;0:0.5:1];
ytl = {'0%','0.3%','0.6%';'0%','40%','80%';'0%','12%','24%';'0%','50%','100%'};
ylm = [0,0,0,0;0.007,0.8,0.25,1];
t_ind = 1:500:length(t);

figure
for i=1:length(ind_species_plot)
    subplot(3,3,(i-1)*3+1)
    plot(t,mean_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0.0],'linewidth',lw);hold on
    plot(t(t_ind),mean_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0.0],'linewidth',lw,'markersize',ms); hold on
    plot(t,mean_x_g2(:,ind_species_plot(i)),'-','color',[0,0.4,0.0],'linewidth',lw);hold on
    plot(t(t_ind),mean_x_g2(t_ind,ind_species_plot(i)),'+','color',[0,0.4,0.0],'linewidth',lw,'markersize',ms);hold on
    plot(t,mean_x_g3(:,ind_species_plot(i)),'-','color',[0,0.1,0.0],'linewidth',lw);hold on
    plot(t(t_ind),mean_x_g3(t_ind,ind_species_plot(i)),'o','color',[0,0.1,0.0],'linewidth',lw,'markersize',ms);hold on
    h1 = plot(t,mean_x_f(:,ind_species_plot(i)),'-','color',[0.5,0,0],'linewidth',lw); hold on
    set(gca,'fontsize',fs)
    xlim([0,44000])
    if i==3
        set(gca,'XTick',0:1e4:4e4)
    else
        set(gca,'XTick',[])
    end
    if i==1
        title('mean','fontsize',fs+1)
    end
    
    subplot(3,3,(i-1)*3+2)
    plot(t,var_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0],'linewidth',lw);hold on
    plot(t(t_ind),var_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0.0],'linewidth',lw,'markersize',ms);
    plot(t,var_x_g2(:,ind_species_plot(i)),'-','color',[0,0.4,0],'linewidth',lw);hold on
    plot(t(t_ind),var_x_g2(t_ind,ind_species_plot(i)),'+','color',[0,0.4,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,var_x_g3(:,ind_species_plot(i)),'-','color',[0,0.1,0],'linewidth',lw);hold on
    plot(t(t_ind),var_x_g3(t_ind,ind_species_plot(i)),'o','color',[0,0.1,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,var_x_f(:,ind_species_plot(i)),'-','color',[0.5,0,0],'linewidth',lw); hold on
    
    set(gca,'fontsize',fs)
    xlim([0,44000])
    if i==3
        set(gca,'XTick',0:1e4:4e4)
    else
        set(gca,'XTick',[])
    end
    if i==1
        title('variance','fontsize',fs+1)
    end
    subplot(3,3,i*3)
    
    plot(t,e_var_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0],'linewidth',lw);hold on
    plot(t(t_ind),e_var_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,e_var_x_g2(:,ind_species_plot(i)),'-','color',[0,0.4,0],'linewidth',lw);hold on
    plot(t(t_ind),e_var_x_g2(t_ind,ind_species_plot(i)),'+','color',[0,0.4,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,e_var_x_g3(:,ind_species_plot(i)),'-','color',[0,0.1,0],'linewidth',lw);hold on
    plot(t(t_ind),e_var_x_g3(t_ind,ind_species_plot(i)),'o','color',[0,0.1,0],'linewidth',lw,'markersize',ms);hold on
    set(gca,'fontsize',fs)
    xlim([0,44000])
    ylim([ylm(1,i),ylm(2,i)])
    set(gca,'YTick',yt(i,:))
    set(gca,'YTickLabel',ytl(i,:))
    if i==3
        set(gca,'XTick',0:1e4:4e4)
    else
        set(gca,'XTick',[])
    end
    if i==1
        title('error in variance','fontsize',fs+1)
    end
    
end

for i=1:3
    subplot(3,3,(i-1)*3+1)
    yl = ylabel(System_MM2_f.state.name(ind_species_plot(i)),'fontsize',fs+1);
    
end

set(gcf,'nextPlot','Add')
h = axes;
hx = xlabel('time','fontsize',fs+1);
set(gca,'visible','off')
set(hx,'visible','on')
