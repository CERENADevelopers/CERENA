clear
close all
clc
%% COMPILATION

[exdir,~,~]=fileparts(which('example_model_1.m'));
% compile the model
cvodewrap('model_example_1','example_model_1_syms',exdir)
% add the model to the path
addpath(genpath([strrep(which('cvodewrap.m'),'cvodewrap.m','') 'models/model_example_1']))

%% SIMULATION

% time vector
t = linspace(0,10,20);
p = [0.5;2;0.5;0.5];
k = [4,8,10,4];

options.sensi = 0;
options.cvode_maxsteps = 1e6;
% load mex into memory
sol = simulate_model_example_1(t,log10(p),k,options);

tic
sol = simulate_model_example_1(t,log10(p),k,options);
disp(['Time elapsed with cvodes: ' num2str(toc) ])

%% ODE15S

ode_system = @(t,x,p,k) [-p(1)*heaviside(t-p(4))*x(1);
    +p(2)*x(1)*exp(-0.1*t)-p(3)*x(2);
    -1.5*x(3)];
% event_fn = @(t,x) [x(3) - x(2);
%     x(3) - x(1)];
% 'Events',event_fn
options_ode15s = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e4);

tic
[~, X_ode15s] = ode15s(@(t,x) ode_system(t,x,p,k),t,k(1:3),options_ode15s);
disp(['Time elapsed with ode15s: ' num2str(toc) ])

%% PLOTTING

figure
c_x = get(gca,'ColorOrder');
subplot(2,2,1)
for ix = 1:size(sol.x,2)
    plot(t,sol.x(:,ix),'.-','Color',c_x(ix,:))
    hold on
    plot(t,X_ode15s(:,ix),'d','Color',c_x(ix,:))
end
stem(sol.root(:,1),sol.root(:,1)*0+10,'r')
stem(sol.root(:,2),sol.root(:,2)*0+10,'k')
legend('x1','x1_{ode15s}','x2','x2_{ode15s}','x3','x3_{ode15s}','x3==x2','x3==x1')
xlabel('time t')
ylabel('x')
box on
subplot(2,2,2)
plot(t,sol.x-X_ode15s,'r--')
legend('error x1','error x2','error x3')

subplot(2,2,3)
plot(t,sol.y,'.-','Color',c_x(1,:))
hold on
plot(t,p(4)*sum(X_ode15s,2),'d','Color',c_x(1,:))
legend('y1','y1_{ode15s}')
xlabel('time t')
ylabel('y')
box on

subplot(2,2,4)
plot(t,sol.y-p(4)*sum(X_ode15s,2),'r--')
legend('error y1')
xlabel('time t')
ylabel('y')
box on


%% FORWARD SENSITIVITY ANALYSIS

options.sensi = 1;

sol = simulate_model_example_1(t,log10(p),k,options);

%% FINITE DIFFERENCES

eps = 1e-4;
xi = log10(p);
for ip = 1:4;
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solp = simulate_model_example_1(t,xip,k,options);
    sx_fd(:,:,ip) = (solp.x - sol.x)/eps;
    sy_fd(:,:,ip) = (solp.y - sol.y)/eps;
    sroot_fd(:,:,ip) = (solp.root - sol.root)/eps;
end




%% PLOTTING
figure
for ip = 1:4
    subplot(4,2,ip*2-1)
    hold on
    for ix = 1:size(sol.x,2)
        plot(t,sol.sx(:,ix,ip),'.-','Color',c_x(ix,:))
        plot(t,sx_fd(:,ix,ip),'d','Color',c_x(ix,:))
    end
    legend('x1','x1_{fd}','x2','x2_{fd}','x3','x3_{fd}')
    title(['state sensitivity for p' num2str(ip)])
    xlabel('time t')
    ylabel('x')
    box on
    
    subplot(4,2,ip*2)
    plot(t,abs(sol.sx(:,:,ip)-sx_fd(:,:,ip)),'r--')
    legend('error x1','error x2','error x3')
    title(['state sensitivity for p' num2str(ip)])
    xlabel('time t')
    ylabel('error')
    set(gca,'YScale','log')
    box on
end

figure
for ip = 1:4
    subplot(4,2,ip*2-1)
    hold on
    for iy = 1:size(sol.y,2)
        plot(t,sol.sy(:,iy,ip),'.-','Color',c_x(iy,:))
        plot(t,sy_fd(:,iy,ip),'d','Color',c_x(iy,:))
    end
    legend('y1','y1_fd')
    title(['observable sensitivity for p' num2str(ip)])
    xlabel('time t')
    ylabel('y')
    box on
    
    subplot(4,2,ip*2)
    plot(t,abs(sol.sy(:,:,ip)-sy_fd(:,:,ip)),'r--')
    legend('error y1')
    title(['error observable sensitivity for p' num2str(ip)])
    xlabel('time t')
    ylabel('error')
    set(gca,'YScale','log')
    box on
end

figure
for ip = 1:4
subplot(4,2,2*ip-1)
bar(1:6,sol.sroot(1:6,:,ip),0.8)
hold on
bar(1:6,sroot_fd(1:6,:,ip),0.4)
legend('x3==x2','x3==x1','x3==x2 fd','x3==x1 fd')
title(['event sensitivity for p' num2str(ip)])
xlabel('event #')
ylabel('y')
box on

subplot(4,2,2*ip)
bar(1:6,sol.sroot(1:6,:,ip)-sroot_fd(1:6,:,ip),0.8)
legend('error x3==x2','error x3==x1')
title(['error event sensitivity for p' num2str(ip)])
xlabel('event #')
ylabel('y')
box on
end



