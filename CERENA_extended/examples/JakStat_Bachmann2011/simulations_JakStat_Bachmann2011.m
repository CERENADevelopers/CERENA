clear
clc
close all
modelDefName = 'modelDef_JakStat_MSB';
%%
load('thetaBestFit.mat') % Loading the parameter values
theta(end) = 1.25e-6;
load('JS_MSB_MM2_f_a_MEC_2_LD_1_a_f_2_System.mat')
System_MM2_f = System;
kappa = zeros(size(System_MM2_f.kappa.variable));
t = 0:0.1:240;
scaleRNA_CIS = 1e4;
tmp = find(strcmp(System_MM2_f.parameter.name,'CISRNATurn'));
theta(tmp) = theta(tmp)/scaleRNA_CIS;
%% Derivation of MA and sMA equations of degree 1, 2 and 3;
%%%% This will take a while; instead you can use the pregenerated models
%%%% below.

% %% full MM simulation
% modelName = 'JS_MSB_MM2_f_a';
% method = 'MEC_2_LD_1_a_f_2';
% System_MM2_f = genmexp(modelName,modelDefName,method);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_f done')
% %% reduced MM with zero closure
% modelName = 'JS_MSB_MM2_g_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g = genmexp(modelName,modelDefName,method);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g done')
% %% reduced MM order 2 with zero closure
% modelName = 'JS_MSB_MM2_g2_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g2.reduction_order = 2;
% System_MM2_g2 = genmexp(modelName,modelDefName,method,System_MM2_g2);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g2 done')
% %% reduced MM order 3 with zero closure
% modelName = 'JS_MSB_MM2_g3_a';
% method = 'MEC_2_LD_1_a_g_2';
% System_MM2_g3.reduction_order = 3;
% System_MM2_g3 = genmexp(modelName,modelDefName,method,System_MM2_g3);
% amiwrap(modelName,[method,'_',modelName,'_syms'])
% disp('MM2_g3 done')
%%
load('JS_MSB_MM2_g_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g = System;
load('JS_MSB_MM2_g2_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g2 = System;
load('JS_MSB_MM2_g3_a_MEC_2_LD_1_a_g_2_System.mat')
System_MM2_g3 = System;

%% simulations
System_MM2_f.sol = simulate_JS_MSB_MM2_f_a(t,theta,kappa);
System_MM2_g.sol = simulate_JS_MSB_MM2_g_a(t,theta,kappa);
System_MM2_g2.sol = simulate_JS_MSB_MM2_g2_a(t,theta,kappa);
System_MM2_g3.sol = simulate_JS_MSB_MM2_g3_a(t,theta,kappa);

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

ind_species_plot= [5,10,18];

yt = [0:0.003:0.006;0:0.15:0.3;0:0.0012:0.0024;0:0.01:0.02];
ytl = {'0%','0.3%','0.6%';'0%','15%','30%';'0%','0.12%','0.24%';'0%','1%','2%'};
ylm = [300,2e4,2500,3000;2500,7e6,1e5,5e5;-0.0002,-0.005,-0.00005,-0.0005;0.006,0.3,0.0025,0.022];
yLabel={'p12EpoR:pJAK2','pSTAT5','CIS'};
t_ind = 1:40:length(t);

figure
for i=1:length(ind_species_plot)
    subplot(3,3,(i-1)*3+1)
    plot(t,mean_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0.0],'linewidth',lw);hold on
    plot(t(t_ind),mean_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0.0],'linewidth',lw,'markersize',ms); hold on
    h1 = plot(t,mean_x_f(:,ind_species_plot(i)),'-','color',[0.5,0,0],'linewidth',lw); hold on
   
    set(gca,'fontsize',fs)
    ylim([-10,ylm(1,i)])
    xlim([0,60])
    if i==3
    set(gca,'XTick',0:20:60)
    else
        set(gca,'XTick',[])
    end
    if i==1
        title('mean','fontsize',fs+1)
    end
    
    subplot(3,3,(i-1)*3+2)
    plot(t,var_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0],'linewidth',lw);hold on
    plot(t(t_ind),var_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0.0],'linewidth',lw,'markersize',ms);
    plot(t,var_x_f(:,ind_species_plot(i)),'-','color',[0.5,0,0],'linewidth',lw); hold on
    
    set(gca,'fontsize',fs)
    xlim([0,60])
    ylim([-30,ylm(2,i)])
  if i==3
    set(gca,'XTick',0:20:60)
    else
        set(gca,'XTick',[])
  end
  if i==1
        title('variance','fontsize',fs+1)
    end
    if i==4
        set(gca,'YTick',0:2.5e5:5e5)
    end
    subplot(3,3,i*3)
    
     plot(t,e_var_x_g(:,ind_species_plot(i)),'-','color',[0,0.7,0],'linewidth',lw);hold on
    plot(t(t_ind),e_var_x_g(t_ind,ind_species_plot(i)),'*','color',[0,0.7,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,e_var_x_g2(:,ind_species_plot(i)),'-','color',[0,0.4,0],'linewidth',lw);hold on
    plot(t(t_ind),e_var_x_g2(t_ind,ind_species_plot(i)),'+','color',[0,0.4,0],'linewidth',lw,'markersize',ms);hold on
    plot(t,e_var_x_g3(:,ind_species_plot(i)),'-','color',[0,0.1,0],'linewidth',lw);hold on
        plot(t(t_ind),e_var_x_g3(t_ind,ind_species_plot(i)),'o','color',[0,0.1,0],'linewidth',lw,'markersize',ms);hold on
        
set(gca,'fontsize',fs)
    xlim([0,60])
    ylim([ylm(3,i),ylm(4,i)])
    %     ylim([0,0.3])
    set(gca,'YTick',yt(i,:))
    set(gca,'YTickLabel',ytl(i,:))
  if i==3
    set(gca,'XTick',0:20:60)
    else
        set(gca,'XTick',[])
  end
  if i==1
        title('error in variance','fontsize',fs+1)
    end
end

for i=1:3
    subplot(3,3,(i-1)*3+1)
    yl = ylabel(yLabel(i),'fontsize',fs+1);
    pyl = get(yl,'Position');
    set(yl,'Position',[-20,pyl(2),pyl(3)])
end

set(gcf,'nextPlot','Add')
h = axes;
hx = xlabel('time','fontsize',fs+1);
set(gca,'visible','off')
set(hx,'visible','on')
set(hx,'Position',[0.5000   -0.05         0])
