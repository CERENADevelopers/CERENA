clear
clc
close all
%%
modelDefName = 'modelDef_linearChain';
t = linspace(0,40,100);
theta = [0.5;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;1;0.2;0.5];
kappa = [];
%% full MM simulation
modelName = 'linearChain_MM2_f_a';
method = 'MEC_2_LD_1_a_f_2';
System_MM2_f = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%% reduced MM with zero closure
modelName = 'linearChain_MM2_g_a';
method = 'MEC_2_LD_1_a_g_2';
System_MM2_g = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%% reduced MM order 3 with zero closure
modelName = 'linearChain_MM2_g2_a';
method = 'MEC_2_LD_1_a_g_2';
System_MM2_g2.reduction_order = 2;
System_MM2_g2 = genmexp(modelName,modelDefName,method,System_MM2_g2);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%% reduced MM order 3 with zero closure
modelName = 'linearChain_MM2_g3_a';
method = 'MEC_2_LD_1_a_g_2';
System_MM2_g3.reduction_order = 3;
System_MM2_g3 = genmexp(modelName,modelDefName,method,System_MM2_g3);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%% reduced MM with zero closure
modelName = 'linearChain_MM3_g_a';
method = 'MEC_3_LD_1_a_g_3';
System_MM3_g = genmexp(modelName,modelDefName,method);
amiwrap(modelName,[method,'_',modelName,'_syms'])
%% simulations
System_MM2_f.sol = simulate_linearChain_MM2_f_a(t,theta,kappa);
System_MM2_g.sol = simulate_linearChain_MM2_g_a(t,theta,kappa);
System_MM2_g2.sol = simulate_linearChain_MM2_g2_a(t,theta,kappa);
System_MM2_g3.sol = simulate_linearChain_MM2_g3_a(t,theta,kappa);
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
%%
% load('20000SSA.mat', 'System_SSA')
% mean_x_ssa = System_SSA.sol.mean_x;
% var_x_ssa = System_SSA.sol.var_x;
%%
fs = 10;
lw = 1;
figure;
for i = 1:System_MM2_g.state.number
    subplot(5,2,i)
    plot(t,mean_x_f(:,i),'-','color',[0,0,0.5],'linewidth',lw); hold on
    plot(t,mean_x_g(:,i),'-','color',[0,0.5,0],'linewidth',lw);hold on
    plot(t,mean_x_g2(:,i),'-','color',[0.2,0.5,0.3],'linewidth',lw);
    plot(t,mean_x_g3(:,i),'-','color',[0.5,0.2,0.3],'linewidth',lw);
%     plot(t,mean_x_ssa(:,i),'--','color',[0,1,0],'linewidth',lw);
ylabel(['mean_X_{',num2str(i),'}'],'fontsize',fs)
set(gca,'fontsize',fs)
if i==1
    leg = legend('full','degree 1','degree 2','degree 3');
    set(leg,'fontsize',fs)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 25 25])
% print('-depsc2','-r1000',['figures/mean']);

figure;
for i = 1:System_MM2_g.state.number
    subplot(5,2,i)
    plot(t,var_x_f(:,i),'-','color',[0,0,0.5],'linewidth',lw); hold on
    plot(t,var_x_g(:,i),'-','color',[0,0.5,0],'linewidth',lw);hold on
    plot(t,var_x_g2(:,i),'-','color',[0.2,0.5,0.3],'linewidth',lw);
    plot(t,var_x_g3(:,i),'-','color',[0.5,0.2,0.3],'linewidth',lw);
%     plot(t,var_x_ssa(:,i),'--','color',[0,1,0],'linewidth',lw);
ylabel(['var_X_{',num2str(i),'}'],'fontsize',fs)
set(gca,'fontsize',fs)
if i==1
    leg = legend('full','degree 1','degree 2','degree 3');
    set(leg,'fontsize',fs)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 25 15])
% print('-depsc2','-r1000',['figures/var']);

%%
figure;
for i = 1:System_MM2_f.state.number
    subplot(5,2,i)
    plot(t,e_mean_x_g(:,i),'-','color',[0,0.5,0],'linewidth',lw); hold on
    plot(t,e_mean_x_g2(:,i),'-','color',[0.7,0.,0.9],'linewidth',lw); hold on
    plot(t,e_mean_x_g3(:,i),'-','color',[0,0.5,0.5],'linewidth',lw); hold on
ylabel(['e_mean_X_{',num2str(i),'}'],'fontsize',fs)
set(gca,'fontsize',fs)
if i==1
    leg = legend('degree 1','degree 2','degree 3');
    set(leg,'fontsize',fs)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 25 15])
% print('-depsc2','-r1000',['figures/e_mean']);

figure;
for i = 1:System_MM2_f.state.number
    subplot(5,2,i)
    plot(t,e_var_x_g(:,i),'-','color',[0,0.5,0],'linewidth',lw); hold on
    plot(t,e_var_x_g2(:,i),'-','color',[0.7,0.,0.9],'linewidth',lw);
    plot(t,e_var_x_g3(:,i),'-','color',[0,0.5,0.5],'linewidth',lw);
ylabel(['e_var_X_{',num2str(i),'}'],'fontsize',fs)
set(gca,'fontsize',fs)
if i==1
    leg = legend('degree 1','degree 2','degree 3');
    set(leg,'fontsize',fs)
end

end
subplot(5,2,9)
xlabel('time','fontsize',fs)
subplot(5,2,10)
xlabel('time','fontsize',fs)
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 25 15])
% print('-depsc2','-r1000',['figures/e_var']);