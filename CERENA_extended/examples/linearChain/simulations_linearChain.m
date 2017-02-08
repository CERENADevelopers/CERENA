clear
clc
close all
modelDefName = 'modelDef_linearChain';
t = linspace(0,40,100);
theta = [0.5;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;0.5];
kappa = [];
%% Simulation using SSA
% Nssa = 20000;
% options.mode = 'constant';
% options.scale = 'absolute';
% System_SSA = simulateSSA_matlab(modelDefName,t,theta,kappa,Nssa,options);
% mean_x_ssa = System_SSA.sol.mean_x;
% var_x_ssa = System_SSA.sol.var_x;
% options.plot_xssa = 0;
% plotSSA(System_SSA,options)
%% Simulation using second-order moment equations (MM2) with low dispersion closure - FULL
modelName = 'linearChain_MM2_f';
method = 'MEC_2_LD_1_a_f';
System_MM2_f = genSimFile(modelName,modelDefName,method);
System_MM2_f.sol = simulate_linearChain_MM2_f(t,theta,kappa);
nSim = 10;
t_tmp = zeros(nSim,1);
for i = 1:nSim
    t1 = tic;
    System_MM2_f.sol = simulate_linearChain_MM2_f(t,theta,kappa);
    t_tmp(i) = toc(t1);
end
tf = mean(t_tmp);
ind_mean_x_f = find(sum(System_MM2_f.MM.sym.state.order>=1,2) == 1);
ind_var_x_f = find((System_MM2_f.MM.sym.state.order(:,end-1)== System_MM2_f.MM.sym.state.order(:,end)).*(sum(System_MM2_f.MM.sym.state.order~=0,2)==2));
mean_x_f = System_MM2_f.sol.x(:,ind_mean_x_f);
var_x_f = System_MM2_f.sol.x(:,ind_var_x_f);
e_mean_x_f = abs(mean_x_f-mean_x_ssa)./mean_x_ssa;
e_var_x_f = abs(var_x_f-var_x_ssa)./var_x_ssa;
%% Simulation using second-order moment equations (MM2) with low dispersion closure - REDUCED
modelName = 'linearChain_MM2_redNorm';
method = 'MEC_2_LD_1_a_n';
System_MM2_n.MM = System_MM2_f.MM;
System_MM2_n = genSimFile(modelName,modelDefName,method,System_MM2_n);
System_MM2_n.sol = simulate_linearChain_MM2_redNorm(t,theta,kappa);
t_tmp = zeros(nSim,1);
for i = 1:nSim
    t1 = tic;
    System_MM2_n.sol = simulate_linearChain_MM2_redNorm(t,theta,kappa);
    t_tmp(i) = toc(t1);
end
tn = mean(t_tmp);
ind_mean_x_n = find(sum(System_MM2_n.MM.sym.state.order>=1,2) == 1);
ind_var_x_n = find((System_MM2_n.MM.sym.state.order(:,end-1)== System_MM2_n.MM.sym.state.order(:,end)).*(sum(System_MM2_n.MM.sym.state.order~=0,2)==2));
mean_x_n = System_MM2_n.sol.x(:,ind_mean_x_n);
var_x_n = System_MM2_n.sol.x(:,ind_var_x_n);
e_mean_x_n = abs(mean_x_n-mean_x_ssa)./mean_x_ssa;
e_var_x_n = abs(var_x_n-var_x_ssa)./var_x_ssa;

e_var_x_n_f = abs(var_x_n-var_x_f)./var_x_f;

%% Simulation using second-order moment equations (MM2) with low dispersion closure - REDUCED, special closure
modelName = 'linearChain_MM2_redSpec';
method = 'MEC_2_LD_1_a_s';
System_MM2_s = genSimFile(modelName,modelDefName,method);
System_MM2_s.sol = simulate_linearChain_MM2_redSpec(t,theta,kappa);
t_tmp = zeros(nSim,1);
for i = 1:nSim
    t1 = tic;
    System_MM2_s.sol = simulate_linearChain_MM2_redSpec(t,theta,kappa);
    t_tmp(i) = toc(t1);
end
ts = mean(t_tmp);
ind_mean_x_s = find(sum(System_MM2_s.MM.sym.state.order>=1,2) == 1);
ind_var_x_s = find((System_MM2_s.MM.sym.state.order(:,end-1)== System_MM2_s.MM.sym.state.order(:,end)).*(sum(System_MM2_s.MM.sym.state.order~=0,2)==2));
mean_x_s = System_MM2_s.sol.x(:,ind_mean_x_s);
var_x_s = System_MM2_s.sol.x(:,ind_var_x_s);
e_mean_x_s = abs(mean_x_s-mean_x_ssa)./mean_x_ssa;
e_var_x_s = abs(var_x_s-var_x_ssa)./var_x_ssa;
%%
lw = 2;
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(5,2,i)
    plot(t,e_mean_x_f(:,i),'-+','color',[0.5,0,0],'linewidth',lw); hold on
    plot(t,e_mean_x_n(:,i),'-o','color',[0,0.5,0],'linewidth',lw);
    plot(t,e_mean_x_s(:,i),'-','color',[0,0,0.5],'linewidth',lw);
end
%%
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(4,3,i)
    plot(t,e_var_x_f(:,i),'-','color',[0.5,0,0],'linewidth',lw); hold on
    plot(t,e_var_x_n(:,i),'-','color',[0,0.5,0],'linewidth',lw);
%     plot(t,e_var_x_s(:,i),'-','color',[0,0,0.5],'linewidth',lw);
ylabel(['e_{var\_X_{',num2str(i),'}}'],'fontsize',fs)
if i==1
    leg = legend('Full','Reduced');
    set(leg,'fontsize',12)
end

end
subplot(4,3,10)
xlabel('time','fontsize',fs)
subplot(4,3,11)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 20])
% print('-depsc2','-r1000',['e_var']);
%%
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(5,2,i)
    plot(t,e_var_x_n_f(:,i),'-','color',[0,0.5,0],'linewidth',lw);
%     plot(t,e_var_x_s(:,i),'-','color',[0,0,0.5],'linewidth',lw);
xlim([0,25])
ylabel(['e_{var\_X_{',num2str(i),'}}'],'fontsize',fs)
if i==1
    leg = legend('Full','Reduced');
    set(leg,'fontsize',12)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 20])
% print('-depsc2','-r1000',['e_var']);
%%
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(5,2,i)
    plot(t,var_x_ssa(:,i),'k+','linewidth',lw); hold on
    plot(t,var_x_f(:,i),'-','color',[0,0,0.5],'linewidth',lw);
    plot(t,var_x_n(:,i),'-','color',[0,0.5,0],'linewidth',lw);
%     plot(t,var_x_s(:,i),'-','color',[0,0,0.5],'linewidth',lw);
ylabel(['X_{',num2str(i),'}'],'fontsize',fs)
if i==1
    leg = legend('SSA','Full','Reduced');
    set(leg,'fontsize',12)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 20])
% print('-depsc2','-r1000',['var']);

%%
lw = 2;
fs = 16;
figure
plot([1,2],[tf,tn],'s-','color',[0.5,0,0],'linewidth',lw)
ylabel('runtime (sec)','fontsize',fs)
xlim([0.5,2.5])
set(gca,'XTick',[1,2])
set(gca,'XTickLabels',{'Full','Reduced'},'fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 7])
print('-depsc2','-r1000',['runtime']);
%%
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(5,2,i)
    plot(t,e_var_x_f(:,i),'-','color',[0.,0,0.5],'linewidth',lw); hold on
    plot(t,e_var_x_n(:,i),'-','color',[0,0.5,0],'linewidth',lw);
    plot(t,e_var_x_s(:,i),'--','color',[0.5,0,0],'linewidth',lw);
ylabel(['e_{var\_X_{',num2str(i),'}}'],'fontsize',fs)
xlim([0,60])
if i==1
    leg = legend('Full','Zero','Independence');
    set(leg,'fontsize',12,'location','northwest')
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 20])
print('-depsc2','-r1000',['e_var_s']);
%%
figure;
for i = 1:size(e_mean_x_f,2)
    subplot(5,2,i)
    plot(t,var_x_ssa(:,i),'k+','linewidth',lw); hold on
    plot(t,var_x_f(:,i),'-','color',[0,0,0.5],'linewidth',lw);
    plot(t,var_x_n(:,i),'-','color',[0,0.5,0],'linewidth',lw);
    plot(t,var_x_s(:,i),'--','color',[0.5,0,0],'linewidth',lw);
ylabel(['X_{',num2str(i),'}'],'fontsize',fs)
xlim([0,60])

if i==1
    leg = legend('SSA','Full','Zero','Independence','location','northwest');
    set(leg,'fontsize',12)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 20])
print('-depsc2','-r1000',['var_s']);

%%
lw = 2;
fs = 16;
figure
plot([1,2,3],[tf,tn,ts],'s-','color',[0.5,0,0],'linewidth',lw)
ylabel('runtime (sec)','fontsize',fs)
xlim([0.5,3.5])
set(gca,'XTick',[1,2,3])
set(gca,'XTickLabels',{'Full','Zero','Independence'},'fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 10 7])
print('-depsc2','-r1000',['runtime_s']);