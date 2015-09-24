function [X,T] = simulateSSAcontinuous(v,S,Ts,x0)

%% SIMULATION
% Initialization
T = Ts(1);
t = Ts(1);
X = x0;
x = x0;
% Simulate till next time point
while t < Ts(end)
    % Evaluation of propensities
    w_i = v(x);
    % w_i = SSA_TF(x); 
    cs_w_i = cumsum(w_i);
    % Sampling of reaction time
    t = 1/cs_w_i(end)*log(1/rand) + t;
    % Sampling of reaction index
    I = max(find([0;cs_w_i(1:end-1)] <= rand*cs_w_i(end)));
    % Update of state
    x = x + S(:,I);
    % Save time point
    X = [X,x];       
    T = [T,t];       
end
